/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNC1DPROPERTYEDITOR_H
#define VRN_TRANSFUNC1DPROPERTYEDITOR_H

#include "voreen/qt/widgets/property/transfunc/transfuncpropertyeditorbase.h"

#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include <QSpinBox>
#include <QLabel>
#include <QComboBox>

namespace voreen {

    class TransFunc1DProperty;

class TransFunc1DPropertyEditor : public TransFuncPropertyEditorBase {
    Q_OBJECT
public:
    /** Constructor.
     * @param prop the transfer function property
     * @param parent the parent frame to which this is added
     */
    TransFunc1DPropertyEditor(TransFunc1DProperty* prop, QWidget* parent = 0);
    /** Destructor */
    ~TransFunc1DPropertyEditor();

    /** @override TransFuncPropertyEditorBase */
    virtual void updateFromProperty();
    //----------------------
    //  Layout functions
    //----------------------
protected:
    /** @override TransFuncpropertyEditorBase */
    virtual QGroupBox* createColorMapSettingsBox();
    /** @override TransFuncpropertyEditorBase */
    virtual QGroupBox* createDomainAndThresholdBox();

    QToolButton* alphaButton_;          ///< Button to change alpha settings
    QMenu* alphaMenu_;                  ///< Menu containing the alpha options
    QDoubleSpinBox* gammaSpin_;         ///< spinbox for gamma value

    QLabel* lowerVolumeBoundLabel_;                 ///< label of lower volume bound
    QLabel* upperVolumeBoundLabel_;                 ///< label of opper volume bound
    QDoubleSpinBox* lowerDomainSpin_;   ///< spinbox for lower mapping value
    QDoubleSpinBox* upperDomainSpin_;   ///< spinbox for upper mapping value
    QDoubleSpinBox* lowerThresholdSpin_;///< spinbox for lower mapping value
    QDoubleSpinBox* upperThresholdSpin_;///< spinbox for upper mapping value
    QSpinBox* cutoffDomainFitting_;     ///< fit domain to data uses min/max value reduced by this percentage
    QComboBox* domainFittingStrategy_;  ///< the fitting stratagy inherited from tfproperty
    QToolButton* fitDomainToData_;      ///< button triggering the domain fitting

    //----------------------
    //  Layout slots
    //----------------------
protected slots:
    /** Applies the button state to the tf */
    void alphaButtonClicked(QAction* action);
    /** Applies the domain from data. */
    void fitDomainClicked();
    /** Applies chnges to the tf. */
    void domainFittingStrategyChanged(int index);
    /** Applies gamma changes to the tf */
    void gammaSpinChanged(double gamma);
    /** Applies domain changes to the tf. */
    void lowerDomainSpinChanged(double value);
    /** Applies domain changes to the tf. */
    void upperDomainSpinChanged(double value);
    /** Applies threshold changes to the tf */
    void lowerThresholdSpinChanged(double value);
    /** Applies threshold changes to the tf */
    void upperThresholdSpinChanged(double value);

    //----------------------
    //  Update functions
    //----------------------
protected:
    /** Updates the alpha button without signals. */
    void updateAlphaButton(TransFuncBase::AlphaMode mode);
    /** Trys to update the volume bound labels. */
    void updateVolumeDataBoundsFromProperty();
    /** Updates the dmoain and threhold values. */
    void updateDomainAndThresholdFromProperty();

    //--------------
    //  Member
    //--------------
protected:
    TransFunc1DProperty* transFunc1DProp_;   ///< pointer to the associated property
    static const std::string loggerCat_; ///< the logger category
};

} // namespace voreen

#endif // VRN_TRANSFUNC1DPROPERTYEDITOR_H
