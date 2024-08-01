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

#ifndef VRN_TRANSFUNCEDITORBASE_H
#define VRN_TRANSFUNCEDITORBASE_H

#include "voreen/qt/widgets/property/transfunc/utils/colorpicker.h"
#include "voreen/qt/widgets/property/transfunc/utils/colorluminancepicker.h"

#include <QWidget>
#include <QToolButton>
#include <QGroupBox>

namespace voreen {

    class TransFuncPropertyBase;

/**
 * Abstract base class for all transfer function widgets. It provides methods to open a Filedialog
 * for loading and saving of a transfer function and a slot for switching coarseness mode on and off.
 */
class TransFuncPropertyEditorBase : public QWidget {
    Q_OBJECT
public:
    /** Constructor */
    TransFuncPropertyEditorBase(TransFuncPropertyBase* property, QWidget* parent = 0);
     /** Destructor */
    virtual ~TransFuncPropertyEditorBase();

    /**
     * The editor must be initialized.
     * This function will call createComponents
     */
    virtual void initialize();

    //----------------------
    //  Layout functions
    //----------------------
protected:
    /** Function to layout all components. */
    void layoutComponents();
    /** Layout the default right components. */
    QWidget* layoutRightComponents();
    /** Must be overridden in sub-class. */
    virtual QWidget* layoutLeftComponents() = 0;


    QGroupBox* createBaseButtonBox();
    /** Override in 2dprimitives */
    virtual QGroupBox* createColorPickerBox();
    /** Override in 1d/2d editor */
    virtual QGroupBox* createColorMapSettingsBox() = 0;
    /** Override in 1d/2d editor */
    virtual QGroupBox* createDomainAndThresholdBox() = 0;
    /** Override if a tooltip widget shall be displayed in the upper right corner */
    virtual QWidget* createTooltip();

    //base buttons
    QToolButton* loadButton_;               ///< button for loading a transfer function
    QToolButton* saveButton_;               ///< button for saving a transfer function
    QToolButton* clearButton_;              ///< button for resetting transfer function to default
    //color picker
    ColorPicker* colorPicker_;              ///< picker for choosing the color of a key in transfer function
    ColorLuminancePicker* colorLumPicker_;  ///< picker for choosing the alpha value of a key


    //----------------------
    //  Layout slots
    //----------------------
protected slots:
    /** Restes the transfer functions. */
    void clearButtonClicked();
    /** Loads a new transfer function. */
    void loadButtonClicked();
    /** Saves the transfer function. */
    void saveButtonClicked();

    //----------------------
    //  General functions
    //----------------------
public:
    /** Used to update the editor to the current proprety state. */
    virtual void updateFromProperty() = 0;

public slots:
    /** Starts or stops the interaction mode. */
    void toggleInteractionMode(bool on);

    //--------------
    //  Member
    //--------------
protected:
    TransFuncPropertyBase* baseProperty_;   ///< pointer to the associated property
    bool initialized_;                      ///< stores, if the editor has been initialized
};

} //namespace voreen

#endif // VRN_TRANSFUNCEDITORBASE_H
