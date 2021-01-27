/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNC1DKEYSPROPERTYEDITOR_H
#define VRN_TRANSFUNC1DKEYSPROPERTYEDITOR_H

#include "voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertyeditor.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

class QToolButton;
class QCheckBox;

namespace voreen {

    class DoubleSlider;
    class TransFunc1DKeysPropertyEditorCanvas;
    class Histogram1D;

    /**
    * TransFuncEditorIntensity is an editor widget for an intensity transfer function.
    * It provides a mapping canvas which displays the keys of the transfer function and the histogram
    * of the associated dataset, a canvas that shows the texture of the transfer function and control
    * elements to adjust the upper and lower thresholds. Furthermore you can change the color of a key
    * with a color and luminancepicker.
    */
    class TransFunc1DKeysPropertyEditor : public TransFunc1DPropertyEditor {
        Q_OBJECT

    public:
        /** Constructor.
        * @param prop the transfer function property
        * @param parent the parent frame to which this is added
        */
        TransFunc1DKeysPropertyEditor(TransFunc1DKeysProperty* prop, QWidget* parent = 0);
        /** Destructor */
        ~TransFunc1DKeysPropertyEditor();

        /** @override TransFuncPropertyEditorBase */
        virtual void updateFromProperty();

        //----------------------
        //  Layout functions
        //----------------------
    protected:
        /** @override TransFuncpropertyEditorBase */
        virtual QWidget* layoutLeftComponents();

    private:
        TransFunc1DKeysPropertyEditorCanvas* mappingCanvas_;    ///< mapping canvas
        DoubleSlider* thresholdSlider_;            ///< 2 slider for adjusting the thresholds

        QCheckBox* showHistogramCB_;    ///< toggles the histogram
        QCheckBox* showTextureCB_;      ///< toggles the background texture
        QToolButton* makeRampButton_;           ///< button to transform a TF into a ramp
        QToolButton* invertButton_;             ///< button to invert the TF

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
        void makeRampButtonClicked();
        /** Inverts the currently loaded TF. */
        void invertMapButtonClicked();

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

        //----------------------
        // Member
        //----------------------
    private:
        TransFunc1DKeysProperty* transFunc1DKeysProp_;
        static const std::string loggerCat_; ///< the logger category
    };

} // namespace voreen

#endif // VRN_TRANSFUNC1DKEYSPROPERTYEDITOR_H

