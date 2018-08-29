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

#ifndef VRN_PLOTPREDICATEDIALOG_H
#define VRN_PLOTPREDICATEDIALOG_H

#include "../datastructures/plotbase.h"
#include "../datastructures/plotcell.h"
#include "../utils/plotpredicatefactory.h"

#include <QDialog>
#include <QComboBox>
#include <QPushButton>
#include <QGridLayout>
#include <QObject>
#include <QLabel>
#include <QDoubleSpinBox>

namespace voreen {

class PlotPredicate;

/**
 * Dialog for editing a PlotPredicate
 **/
class PlotPredicateDialog : public QDialog {
    Q_OBJECT

public:
    /**
     * Constructs a new PlotPredicateDialog.
     *
     * \param   predicate               initial PlotPredicate of the dialog, may be NULL if no Predicate to preselect
     * \param   tagsOnly                if false only numeric predicates and parameters, if true only tag predicates and
     *                                  parameters in this dialog
     * \param   simplePredicatesOnly    if true, only simple predicates are selectable
     * \param   parent                  parent widget
     *
     * \note    This dialog makes a deep copy of \a predicate, so it is safe to delete it at any time after this method returns.
     **/
    PlotPredicateDialog(const PlotPredicate* predicate, bool tagsOnly, bool simplePredicatesOnly = true, QWidget* parent = 0);

    /**
     * Destructor
     **/
    ~PlotPredicateDialog();

    /**
     * Returns the currently selected PlotPredicate.
     *
     * \note    The returned PlotPredicate will be deleted on dialog destruction so make sure to perform
     *          deep copies of it.
     **/
    const PlotPredicate* getPlotPredicate() const;

private:
    /// creates all widgets (includes calling updateInnerWidgets())
    void createWidgets();

    /// updates all inner widgets for setting up selected PlotPredicate
    void updateInnerWidgets();

    PlotPredicate* predicate_;              ///< current PlotPredicate
    bool tagsOnly_;                         ///< flag whether to display only predicates for tags or not
    bool simplePredicatesOnly_;             ///< flag whether to display only simple predicates

    QGridLayout* mainLayout_;               ///< main layout of the dialog
    QComboBox* cbPredicates_;               ///< ComboBox with the predicates

    std::vector<QLabel*> labelVector_;      ///< vector of the labels for the edit fields
    std::vector<QLineEdit*> editVector_;    ///< vector of the edit fields
    std::vector<QDoubleSpinBox*> sbVector_; ///< vector of the SpinBoxes

    QPushButton* btnOK_;                    ///< OK button
    QPushButton* btnCancel_;                ///< Cancel button

private slots:
    void clickedOk();
    void clickedCancel();
    void choiceChange(int index);
    void updatePredicate();

};

}
#endif // VRN_PLOTPREDICATEDIALOG_H
