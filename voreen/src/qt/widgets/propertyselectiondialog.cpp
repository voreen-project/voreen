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

#include "voreen/qt/widgets/propertyselectiondialog.h"

#include "voreen/core/network/processornetwork.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/property.h"

#include <QDialogButtonBox>
#include <QTreeWidget>
#include <QVariant>
#include <QVBoxLayout>

namespace voreen {

namespace {

inline QStringList stringListFromStdString(const std::string& string) {
    return QStringList(QString::fromStdString(string));
}

// reinterpret_cast to store the pointers as integer values in variants
inline QVariant processorToVariant(Processor* processor) {
    return qVariantFromValue((void*)processor);
}

inline Processor* variantToProcessor(const QVariant& variant) {
    return (Processor*)variant.value<void*>();
}

inline QVariant propertyToVariant(Property* property) {
    return qVariantFromValue((void*)property);
}

inline Property* variantToProperty(const QVariant& variant) {
    return (Property*)variant.value<void*>();
}

}

PropertySelectionDialog::PropertySelectionDialog(QList<Property*>& resultList, const ProcessorNetwork* network,
                                                QWidget* parent)
    : QDialog(parent)
    , resultList_(resultList)
    , network_(network)
{
    QVBoxLayout* layout = new QVBoxLayout(this);
    treeWidget_ = new QTreeWidget;
    treeWidget_->setHeaderHidden(true);
    treeWidget_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    layout->addWidget(treeWidget_);


    QDialogButtonBox* buttonBox = new QDialogButtonBox;
    buttonBox->addButton(QDialogButtonBox::Ok);
    buttonBox->addButton(QDialogButtonBox::Close);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
    layout->addWidget(buttonBox, 0, Qt::AlignRight);

    foreach (Processor* processor, network->getProcessors()) {
        if (processor->getProperties().size() > 0) {
            QTreeWidgetItem* processorItem = new QTreeWidgetItem(stringListFromStdString(processor->getGuiName()));
            processorItem->setData(0, Qt::UserRole, processorToVariant(processor));
            foreach (Property* property, processor->getProperties()) {
                QTreeWidgetItem* propertyItem = new QTreeWidgetItem(processorItem, stringListFromStdString(property->getGuiName()));
                propertyItem->setData(0, Qt::UserRole, propertyToVariant(property));
            }
            treeWidget_->addTopLevelItem(processorItem);
        }
    }

    treeWidget_->selectAll();
}

void PropertySelectionDialog::accept() {
    QList<QTreeWidgetItem*> items = treeWidget_->selectedItems();
    foreach (QTreeWidgetItem* item, items) {
        if (item->childCount() == 0) {
            Property* property = variantToProperty(item->data(0, Qt::UserRole));
            resultList_.append(property);
        }
        else {
            QList<QTreeWidgetItem*> children = item->takeChildren();
            foreach (QTreeWidgetItem* childItem, children) {
                Property* property = variantToProperty(childItem->data(0, Qt::UserRole));
                resultList_.append(property);
            }
        }
    }
    QDialog::accept();
}

} // namespace
