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

#include "interactivelistpropertywidget.h"

#include "../../properties/interactivelistproperty.h"

#include "voreen/qt/widgets/customlabel.h"

#include <QApplication>
#include <QListWidget>
#include <QDrag>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QMenu>
#include <QMimeData>

namespace voreen {

class ItemListWidget : public QListWidget {
public:

    ItemListWidget(InteractiveListProperty* property, QWidget* parent)
        : QListWidget(parent)
        , property_(property)
    {
        setAcceptDrops(true);
        setMinimumHeight(150);
    }

    void rebuild() {
        clear();

        const std::vector<std::string>& items = property_->getItems();
        for (int idx : property_->getInputIndices()) {
            QListWidgetItem* item = new QListWidgetItem(QString::fromStdString(items[idx]));
            addItem(item);
        }
    }

    void mousePressEvent(QMouseEvent* event) {
        if (event->button() == Qt::LeftButton) {
            dragStartPosition_ = event->pos();
            QListWidgetItem* selectedItem = itemAt(dragStartPosition_);
            if(selectedItem) {
                setCurrentRow(row(selectedItem));
            }
        }
    }

    void mouseMoveEvent(QMouseEvent* event) {
        if (!(event->buttons() & Qt::LeftButton))
            return;
        if ((event->pos() - dragStartPosition_).manhattanLength()
            < QApplication::startDragDistance())
            return;

        QListWidgetItem* selectedItem = itemAt(dragStartPosition_);
        if(!selectedItem)
            return;

        int itemIdx = property_->getInputIndices()[row(selectedItem)];
        QByteArray data;
        QDataStream dataStream(&data, QIODevice::WriteOnly);
        dataStream << itemIdx;

        QDrag* drag = new QDrag(this);
        QMimeData* mimeData = new QMimeData;
        mimeData->setData("filter/item", data);
        drag->setMimeData(mimeData);
        drag->exec();
    }

    void dragEnterEvent(QDragEnterEvent* event) {
        if (event->mimeData()->hasFormat("filter/instance")) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dragMoveEvent(QDragMoveEvent* event) {
        if (event->mimeData()->hasFormat("filter/instance")) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dropEvent(QDropEvent* event) {
        if(event->mimeData()->hasFormat("filter/instance")) {
            int instanceIdx = 0;
            QByteArray data = event->mimeData()->data("filter/instance");
            QDataStream dataStream(&data, QIODevice::ReadOnly);
            dataStream >> instanceIdx;
            event->acceptProposedAction();
            property_->removeInstance(property_->getInstances()[instanceIdx].getInstanceId());
        }
        else {
            event->ignore();
        }
    }

    void keyPressEvent(QKeyEvent* event) {
        QListWidgetItem* selectedItem = item(currentRow());
        if(event->key() == Qt::Key_Return && selectedItem) {
            property_->addInstance(selectedItem->text().toStdString());
        }
        QListWidget::keyPressEvent(event);
    }

private:

    InteractiveListProperty* property_;
    QPoint dragStartPosition_;
};

class InstanceListWidget : public QListWidget {
public:
    InstanceListWidget(InteractiveListProperty* property, QWidget* parent)
        : QListWidget(parent)
        , property_(property)
    {
        setAcceptDrops(true);
        setMinimumHeight(150);
        setContextMenuPolicy(Qt::CustomContextMenu);
        connect(this, &InstanceListWidget::customContextMenuRequested, this, &InstanceListWidget::showContextMenu);
    }

    void rebuild() {
        clear();

        for (const InteractiveListProperty::Instance& instance : property_->getInstances()) {
            QListWidgetItem* item = new QListWidgetItem(QString::fromStdString(instance.getName()));
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            addItem(item);

            // Strike out item, if not active.
            if(!instance.isActive()) {
                QFont font = item->font();
                font.setStrikeOut(true);
                item->setFont(font);
            }
        }
        setCurrentRow(property_->getSelectedInstance());
    }

    void mousePressEvent(QMouseEvent* event) {
        if (event->button() & (Qt::LeftButton | Qt::RightButton)) {
            dragStartPosition_ = event->pos();
            QListWidgetItem* selectedItem = itemAt(dragStartPosition_);
            if(selectedItem) {
                property_->setSelectedInstance(row(selectedItem));
            }
            else {
                property_->setSelectedInstance(-1);
            }
        }
    }

    void mouseMoveEvent(QMouseEvent* event) {
        if (!(event->buttons() & Qt::LeftButton))
            return;
        if ((event->pos() - dragStartPosition_).manhattanLength()
            < QApplication::startDragDistance())
            return;

        QListWidgetItem* selectedItem = itemAt(dragStartPosition_);
        if(!selectedItem)
            return;

        QByteArray data;
        QDataStream dataStream(&data, QIODevice::WriteOnly);
        dataStream << row(selectedItem);

        QDrag* drag = new QDrag(this);
        QMimeData* mimeData = new QMimeData;
        mimeData->setData("filter/instance", data);
        drag->setMimeData(mimeData);
        drag->exec();
    }

    void dragEnterEvent(QDragEnterEvent* event) {
        if (event->mimeData()->hasFormat("filter/item") || event->mimeData()->hasFormat("filter/instance")) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dragMoveEvent(QDragMoveEvent* event) {
        if (event->mimeData()->hasFormat("filter/item") || event->mimeData()->hasFormat("filter/instance")) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dropEvent(QDropEvent* event) {
        if(event->mimeData()->hasFormat("filter/instance")) {
            int position = row(itemAt(event->pos()));
            if(position == -1) {
                event->ignore();
                return;
            }

            int instanceIdx = 0;
            QByteArray data = event->mimeData()->data("filter/instance");
            QDataStream dataStream(&data, QIODevice::ReadOnly);
            dataStream >> instanceIdx;

            property_->swapInstances(property_->getInstances()[instanceIdx].getInstanceId(), position);
            event->acceptProposedAction();
        }
        else if(event->mimeData()->hasFormat("filter/item")) {
            int itemIdx = 0;
            QByteArray data = event->mimeData()->data("filter/item");
            QDataStream dataStream(&data, QIODevice::ReadOnly);
            dataStream >> itemIdx;

            property_->addInstance(property_->getItems()[itemIdx]);
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void keyPressEvent(QKeyEvent* event) {
        tgtAssert(currentRow() == property_->getSelectedInstance(), "Selection does not match");
        tgtAssert(count() == static_cast<int>(property_->getInstances().size()), "Instance count does not match");

        // Accept beforehand and ignore if necessary.
        event->accept();

        int current = property_->getSelectedInstance();
        if(current == -1 || property_->getInstances().empty()) {
            return;
        }

        InteractiveListProperty::Instance instance = property_->getInstances()[current];
        if(event->key() == Qt::Key_Up) {
            current--;
        }
        else if(event->key() == Qt::Key_Down) {
            current++;
        }
        else if(event->key() == Qt::Key_Delete) {
            property_->removeInstance(instance.getInstanceId());
            return;
        }
        else {
            event->ignore();
            return;
        }

        current = tgt::clamp(current, 0, count()-1);

        if(event->modifiers() & Qt::Modifier::CTRL) {
            property_->swapInstances(instance.getInstanceId(), current);
        }
        else {
            property_->setSelectedInstance(current);
        }
    }

public slots:

    /**
     * This function will show the context menu for the item below the cursor.
     */
    void showContextMenu(const QPoint& pos) {
        QListWidgetItem* selectedItem = itemAt(pos);
        if(selectedItem) {
            QMenu menu;
            // Note: using a lambda function would be good here, but it's supported since Qt version 5.6, only.
            // Ubuntu 16.04 only provides 5.5.1.
            menu.addAction(tr("Rename"), this, &InstanceListWidget::renameItem);
            menu.exec(mapToGlobal(pos));
        }
    }

    /**
     * So this slot has to be created in order for the rename action in the context menu (see above) to be executed.
     */
    void renameItem() {
        QListWidgetItem* selectedItem = itemAt(dragStartPosition_);
        tgtAssert(selectedItem, "item null");
        editItem(selectedItem);
    }

    /**
     * This slot is triggered, when the item got its name, be it after creation or change.
     */
    void textChanged(QListWidgetItem* item) {
        int instanceIdx = row(item);
        InteractiveListProperty::Instance& instance = property_->getInstances()[instanceIdx];

        // Update the instance's name accordingly.
        if(!item->text().isEmpty()) {
            instance.setName(item->text().toStdString());

            // TODO: The state is well defined here, since both instance and widget are up to date.
            //  However, owning processors have to be informed in order to be up to date as well.
            //property_->invalidate();
        }
        // Otherwise, discard.
        else {
            item->setText(QString::fromStdString(instance.getName()));
        }
    }

private:

    InteractiveListProperty* property_;
    QPoint dragStartPosition_;
};


InteractiveListPropertyWidget::InteractiveListPropertyWidget(InteractiveListProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , itemTable_(nullptr)
    , instanceTable_(nullptr)
    , property_(prop)
{
    tgtAssert(property_, "null pointer");

    QBoxLayout* itemLayout = new QVBoxLayout();
    if(!property_->getItemLabel().empty()) {
        QLabel* label = new QLabel(QString::fromStdString(property_->getItemLabel()));
        itemLayout->addWidget(label);
    }

    QBoxLayout* instanceLayout = new QVBoxLayout();
    if(!property_->getInstanceLabel().empty()) {
        QLabel* label = new QLabel(QString::fromStdString(property_->getInstanceLabel()));
        instanceLayout->addWidget(label);
    }

    itemTable_ = new ItemListWidget(property_, nullptr);
    instanceTable_ = new InstanceListWidget(property_, nullptr);

    itemTable_->setToolTip(tr("Use Drag and Drop or press RETURN to add an instance of the selected item"));
    instanceTable_->setToolTip(tr("Use Drag and Drop or CTRL + Up/Down to move items"));

    // Add widgets.
    itemLayout->addWidget(itemTable_);
    instanceLayout->addWidget(instanceTable_);

    addLayout(itemLayout);
    addLayout(instanceLayout);

    // Fill tables.
    updateTable();
}

InteractiveListPropertyWidget::~InteractiveListPropertyWidget() {
}

//----------------------------------------------------------
//  Overrides
//----------------------------------------------------------
void InteractiveListPropertyWidget::updateFromPropertySlot() {
    updateTable();
}

//----------------------------------------------------------
//  Helpers
//----------------------------------------------------------
void InteractiveListPropertyWidget::updateTable() {
    itemTable_->rebuild();
    instanceTable_->rebuild();
}

} // namespace
