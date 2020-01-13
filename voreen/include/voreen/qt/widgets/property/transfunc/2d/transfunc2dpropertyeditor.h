/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNC2DPROPERTYEDITOR_H
#define VRN_TRANSFUNC2DPROPERTYEDITOR_H

#include "voreen/qt/widgets/property/transfunc/transfuncpropertyeditorbase.h"

#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include <QSpinBox>
#include <QLabel>
#include <QComboBox>

namespace voreen {

    class TransFunc2DProperty;

class TransFunc2DPropertyEditor : public TransFuncPropertyEditorBase {
    Q_OBJECT
public:
    /** Constructor.
     * @param prop the transfer function property
     * @param parent the parent frame to which this is added
     */
    TransFunc2DPropertyEditor(TransFunc2DProperty* prop, QWidget* parent = 0);
    /** Destructor */
    ~TransFunc2DPropertyEditor();

    /** @override TransFuncPropertyEditorBase */
    virtual void updateFromProperty();
    /** Should return the gui name of channel x. */
    virtual std::string getXGuiName() = 0;
    /** Should return the gui name of channel y. */
    virtual std::string getYGuiName() = 0;
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
    QDoubleSpinBox* gammaSpinX_;        ///< spinbox for gamma value in X dim
    QDoubleSpinBox* gammaSpinY_;        ///< spinbox for gamma value in Y dim

    QLabel* lowerVolumeBoundLabelX_;    ///< label of lower volume bound in X dim
    QLabel* upperVolumeBoundLabelX_;    ///< label of opper volume bound in Y dim
    QLabel* lowerVolumeBoundLabelY_;    ///< label of lower volume bound in X dim
    QLabel* upperVolumeBoundLabelY_;    ///< label of opper volume bound in Y dim

    QDoubleSpinBox* lowerDomainSpinX_;   ///< spinbox for lower mapping value in X dim
    QDoubleSpinBox* upperDomainSpinX_;   ///< spinbox for upper mapping value in X dim
    QDoubleSpinBox* lowerDomainSpinY_;   ///< spinbox for lower mapping value in Y dim
    QDoubleSpinBox* upperDomainSpinY_;   ///< spinbox for upper mapping value in Y dim
    QDoubleSpinBox* lowerThresholdSpinX_;///< spinbox for lower mapping value in X dim
    QDoubleSpinBox* upperThresholdSpinX_;///< spinbox for upper mapping value in X dim
    QDoubleSpinBox* lowerThresholdSpinY_;///< spinbox for lower mapping value in Y dim
    QDoubleSpinBox* upperThresholdSpinY_;///< spinbox for upper mapping value in Y dim

    QSpinBox* cutoffDomainFittingX_;    ///< fit domain to data uses min/max value reduced by this percentage in X dim
    QSpinBox* cutoffDomainFittingY_;    ///< fit domain to data uses min/max value reduced by this percentage in Y dim
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
    void gammaSpinChanged();
    /** Applies domain changes to the tf. */
    void lowerDomainSpinChanged();
    /** Applies domain changes to the tf. */
    void upperDomainSpinChanged();
    /** Applies threshold changes to the tf */
    void lowerThresholdSpinChanged();
    /** Applies threshold changes to the tf */
    void upperThresholdSpinChanged();

    //----------------------
    //  Update functions
    //----------------------
protected:
    /** Updates the alpha button without signals. */
    void updateAlphaButton(TransFuncBase::AlphaMode mode);
    /**
     * Trys to update the volume bound labels.
     * Pure virtual, since X and Y are not defined in 2d functions.
     */
    virtual void updateVolumeDataBoundsFromProperty() = 0;
    /** Updates the dmoain and threhold values. */
    void updateDomainAndThresholdFromProperty();

    //--------------
    //  Member
    //--------------
protected:
    TransFunc2DProperty* transFunc2DProp_;   ///< pointer to the associated property
    static const std::string loggerCat_; ///< the logger category
};


} // namespace voreen

#endif // VRN_TRANSFUNC2DPROPERTYEDITOR_H
