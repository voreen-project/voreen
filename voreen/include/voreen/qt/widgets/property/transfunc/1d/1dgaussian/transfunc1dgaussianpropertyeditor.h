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

#ifndef VRN_TRANSFUNC1DGAUSSIANPROPERTYEDITOR_H
#define VRN_TRANSFUNC1DGAUSSIANPROPERTYEDITOR_H

#include "voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertyeditor.h"
#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"

class QToolButton;
class QCheckBox;

namespace voreen {

    class DoubleSlider;
    class TransFunc1DGaussianPropertyEditorCanvas;
    class Histogram1D;

    /**
    * TransFunc1DGaussianPropertyEditor is an editor widget for an intensity transfer function.
    * It provides a mapping canvas which displays the curves of the transfer function and the histogram
    * of the associated dataset, a canvas that shows the texture of the transfer function and control
    * elements to adjust the upper and lower thresholds. Furthermore you can change the color of a curve
    * with a color and luminancepicker.
    */
    class TransFunc1DGaussianPropertyEditor : public TransFunc1DPropertyEditor {
        Q_OBJECT

    public:
        /** Constructor.
        * @param prop the transfer function property
        * @param parent the parent frame to which this is added
        */
        TransFunc1DGaussianPropertyEditor(TransFunc1DGaussianProperty* prop, QWidget* parent = 0);
        /** Destructor */
        ~TransFunc1DGaussianPropertyEditor();

        /** @override TransFuncPropertyEditorBase */
        virtual void updateFromProperty();

        //----------------------
        //  Layout functions
        //----------------------
    protected:
        /** @override TransFuncpropertyEditorBase */
        virtual QWidget* layoutLeftComponents();
        virtual QWidget* createTooltip();

    private:
        TransFunc1DGaussianPropertyEditorCanvas* mappingCanvas_;    ///< mapping canvas
        DoubleSlider* thresholdSlider_;            ///< 2 slider for adjusting the thresholds

        QCheckBox* showHistogramCB_;    ///< toggles the histogram
        QCheckBox* showTextureCB_;      ///< toggles the background texture
        QToolButton* autoTransFuncButton;   ///< button to transform a TF into a ramp
        QToolButton* tooltipButton_;    ///< button to show a tooltip

        QCheckBox* computeHistogram_;   ///< auto computation of histogram?

        //----------------------
        //  Layout slots
        //----------------------
        protected slots:
        /** Applies threshold changes to the tf */
        void sliderThresholdChanged();
        /** Toggles the histogram in the mapping canvas. */
        void showTextureToggled(int state);
        /** Toggles the texture in the mapping canvas. */
        void showHistogramToggled(int state);
        /** Transform the currently loaded TF into a ramp. */
        void autoTransFuncClicked();

        void computeHistogramToggled(int state);

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
        *
        */
        void markerColorChanged(int h, int s, int v);

        /**
        * Shows a tooltip if the tooltip button was clicked.
        */
        void showTooltip();

        //----------------------
        // Member
        //----------------------
    private:
        TransFunc1DGaussianProperty* transFunc1DGaussianProp_;
        static const std::string loggerCat_; ///< the logger category
    };

} // namespace voreen

#endif // VRN_TRANSFUNC1DGAUSSIANPROPERTYEDITOR_H

