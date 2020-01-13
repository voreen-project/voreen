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

#include "plotpredicatedialog.h"

#include "../datastructures/plotpredicate.h"
#include "../utils/plotpredicatefactory.h"

#include <QDialogButtonBox>
#include <QLineEdit>
#include <QAction>
#include <limits>

namespace voreen {

PlotPredicateDialog::PlotPredicateDialog(const PlotPredicate* predicate, bool tagsOnly, bool simplePredicatesOnly, QWidget* parent)
    : QDialog(parent)
    , predicate_(0)
    , tagsOnly_(tagsOnly)
    , simplePredicatesOnly_(simplePredicatesOnly)
    , mainLayout_(0)
    , cbPredicates_(0)
    , btnOK_(0)
    , btnCancel_(0)
{
    setWindowTitle(QString(tr("Select Predicates",0)));
    setModal(true);

    if (predicate != 0)
        predicate_ = predicate->clone();

    createWidgets();
}

PlotPredicateDialog::~PlotPredicateDialog() {
    delete predicate_;
}

void PlotPredicateDialog::createWidgets() {
    mainLayout_ = new QGridLayout();
    QLabel* lblPredicates = new QLabel(tr("Select a Predicate:"));
    mainLayout_->addWidget(lblPredicates, 0, 0);

    std::vector<std::string> typeStrings;
    if (tagsOnly_ && simplePredicatesOnly_)
        typeStrings = PlotPredicateFactory::getAllTypeStrings(PlotPredicateFactory::SIMPLE_TAG_PREDICATES_ONLY);
    else if (tagsOnly_ && !simplePredicatesOnly_)
        typeStrings = PlotPredicateFactory::getAllTypeStrings(PlotPredicateFactory::ALL_TAG_PREDICATES_ONLY);
    else if (!tagsOnly_ && simplePredicatesOnly_)
        typeStrings = PlotPredicateFactory::getAllTypeStrings(PlotPredicateFactory::SIMPLE_VALUE_PREDICATES_ONLY);
    else
        typeStrings = PlotPredicateFactory::getAllTypeStrings(PlotPredicateFactory::ALL_VALUE_PREDICATES_ONLY);
    cbPredicates_ = new QComboBox();

    // fill combobox
    int i = 1;
    std::string typeString = (predicate_ != 0) ? PlotPredicateFactory::getInstance()->getSerializableTypeString(typeid(*predicate_)) : "";
    cbPredicates_->addItem(tr("No Predicate"), QVariant(tr("")));
    cbPredicates_->setCurrentIndex(0);
    for (std::vector<std::string>::const_iterator it = typeStrings.begin(); it != typeStrings.end(); ++it, ++i) {
        cbPredicates_->addItem(QString::fromStdString(*it), QVariant(QString::fromStdString(*it)));
        // select item if it matches current predicate
        if (*it == typeString)
            cbPredicates_->setCurrentIndex(i);
    }

    mainLayout_->addWidget(cbPredicates_, 0, 1, 1, 2);

    btnOK_ = new QPushButton(tr("OK"));
    btnCancel_ = new QPushButton(tr("Cancel"));
    updateInnerWidgets();

    setLayout(mainLayout_);

    connect(cbPredicates_, SIGNAL(currentIndexChanged(int)), this, SLOT(choiceChange(int)));
    connect(btnOK_, SIGNAL(clicked()), this, SLOT(clickedOk()));
    connect(btnCancel_, SIGNAL(clicked()), this, SLOT(clickedCancel()));
}

const PlotPredicate* PlotPredicateDialog::getPlotPredicate() const {
    return predicate_;
}

void PlotPredicateDialog::updateInnerWidgets() {
    // remove all widgets first
    for (std::vector<QLabel*>::iterator it = labelVector_.begin(); it != labelVector_.end(); ++it)
        delete *it;
    labelVector_.clear();

    for (std::vector<QDoubleSpinBox*>::iterator it = sbVector_.begin(); it != sbVector_.end(); ++it)
        delete *it;
    sbVector_.clear();

    for (std::vector<QLineEdit*>::iterator it = editVector_.begin(); it != editVector_.end(); ++it)
        delete *it;
    editVector_.clear();

    // now create the new widgets
    int innerWidgetCount = (predicate_ == 0 ? 0 : predicate_->getNumberOfThresholdValues());
    if (innerWidgetCount > 0) {
        std::vector<std::string> titles = predicate_->getThresholdTitles();
        std::vector<PlotCellValue> values = predicate_->getThresholdValues();

        for (int i = 0; i < innerWidgetCount; ++i) {
            labelVector_.push_back(new QLabel(QString::fromStdString(titles[i])));
            mainLayout_->addWidget(labelVector_[i], i+2, 0);

            if (tagsOnly_) {
                editVector_.push_back(new QLineEdit(QString::fromStdString(values[i].getTag())));
                connect(editVector_[i], SIGNAL(editingFinished()), this, SLOT(updatePredicate()));
                mainLayout_->addWidget(editVector_[i], i+2, 1, 1, 2);
            }
            else {
                sbVector_.push_back(new QDoubleSpinBox());
                // REMARK: std::numeric_limits<double>::min() does not work here :( hopefully -1000000000 is low enough
                sbVector_[i]->setMinimum(-1000000000);
                sbVector_[i]->setMaximum(std::numeric_limits<double>::max());
                sbVector_[i]->setDecimals(5);
                sbVector_[i]->setValue(values[i].getValue());
                connect(sbVector_[i], SIGNAL(editingFinished()), this, SLOT(updatePredicate()));
                mainLayout_->addWidget(sbVector_[i], i+2, 1, 1, 2);
            }
        }
    }

    // move dialog buttons to the appropriate position
    mainLayout_->addWidget(btnOK_, innerWidgetCount + 3, 1);
    mainLayout_->addWidget(btnCancel_, innerWidgetCount + 3, 2);
}

void PlotPredicateDialog::choiceChange(int /*index*/) {
    PlotPredicate* copy = predicate_;

    // a dynamic_cast to 0 is also OK here as it says 'no predicate selected'
    predicate_ = dynamic_cast<PlotPredicate*>(PlotPredicateFactory::getInstance()->createSerializableType(cbPredicates_->currentText().toStdString()));

    if (predicate_ != 0) {
        int newThresholdCount = predicate_->getNumberOfThresholdValues();
        std::vector<PlotCellValue> values;
        if (copy != 0) {
            values = copy->getThresholdValues();
            delete copy;
            copy = 0;
        }
        values.resize(newThresholdCount, (tagsOnly_ ? PlotCellValue("") : PlotCellValue(0)));
        predicate_->setThresholdValues(values);

    }
    updateInnerWidgets();
}

void PlotPredicateDialog::updatePredicate() {
    if (predicate_ != 0) {
        int newThresholdCount = predicate_->getNumberOfThresholdValues();
        std::vector<PlotCellValue> values;
        for (int i = 0; i < newThresholdCount; ++i) {
            if (tagsOnly_) {
                values.push_back(PlotCellValue(editVector_[i]->text().toStdString()));
            }
            else {
                values.push_back(PlotCellValue(sbVector_[i]->value()));
            }
        }
        predicate_->setThresholdValues(values);
    }
}

void PlotPredicateDialog::clickedOk() {
    accept();
}

void PlotPredicateDialog::clickedCancel() {
    reject();
}


} // namespace voreen
