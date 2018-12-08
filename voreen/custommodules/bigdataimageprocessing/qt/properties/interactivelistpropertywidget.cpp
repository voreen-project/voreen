/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "interactivelistpropertywidget.h"

#include "../../properties/interactivelistproperty.h"

#include "voreen/qt/widgets/customlabel.h"

#include <QApplication>
#include <QListWidget>
#include <QDragLeaveEvent>
#include <QDropEvent>
#include <QMimeData>
#include <QDrag>
#include <QAbstractItemModel>

namespace voreen {

class ItemListWidget : public QListWidget {
public:

    ItemListWidget(InteractiveListProperty* property, QWidget* parent)
        : QListWidget(parent)
        , property_(property)
        , partner_(nullptr)
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

        QDrag* drag = new QDrag(this);
        QMimeData* mimeData = new QMimeData;
        mimeData->setText(selectedItem->text());
        drag->setMimeData(mimeData);
        drag->exec();
    }

    void dragEnterEvent(QDragEnterEvent* event) {
        if (event->mimeData()->hasText() && event->source() == partner_) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dragMoveEvent(QDragMoveEvent* event) {
        if (event->mimeData()->hasText() && event->source() == partner_ ) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dropEvent(QDropEvent* event) {
        if(event->mimeData()->hasText() && event->source() == partner_) {
            //event->setDropAction(Qt::CopyAction);
            event->acceptProposedAction();
            property_->removeInstance(event->mimeData()->text().toStdString());
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

    void setPartner(QWidget* partner) {
        partner_ = partner;
    }

private:

    InteractiveListProperty* property_;
    QWidget* partner_;
    QPoint dragStartPosition_;
};

class InstanceListWidget : public QListWidget {
public:
    InstanceListWidget(InteractiveListProperty* property, QWidget* parent)
        : QListWidget(parent)
        , property_(property)
        , partner_(nullptr)
    {
        setAcceptDrops(true);
        setMinimumHeight(150);
    }

    void rebuild() {
        clear();

        for (const InteractiveListProperty::Instance& instance : property_->getInstances()) {
            QListWidgetItem* item = new QListWidgetItem(QString::fromStdString(instance.name_));
            addItem(item);
        }
        setCurrentRow(property_->getSelectedInstance());
    }

    void mousePressEvent(QMouseEvent* event) {
        if (event->button() == Qt::LeftButton) {
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

        QDrag* drag = new QDrag(this);
        QMimeData* mimeData = new QMimeData;
        mimeData->setText(selectedItem->text());
        drag->setMimeData(mimeData);
        drag->exec();
    }

    void dragEnterEvent(QDragEnterEvent* event) {
        if (event->mimeData()->hasText() && (event->source() == partner_ || event->source() == this)) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dragMoveEvent(QDragMoveEvent* event) {
        if (event->mimeData()->hasText() && (event->source() == partner_ || event->source() == this)) {
            event->acceptProposedAction();
        }
        else {
            event->ignore();
        }
    }

    void dropEvent(QDropEvent* event) {
        if(event->mimeData()->hasText()) {
            if(event->source() == this) {
                int position = row(itemAt(event->pos()));
                if(position == -1) {
                    event->ignore();
                    return;
                }
                property_->moveInstance(event->mimeData()->text().toStdString(), position);
            }
            else {
                property_->addInstance(event->mimeData()->text().toStdString());
            }

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
            property_->removeInstance(instance.instanceId_);
            return;
        }
        else {
            event->ignore();
            return;
        }

        current = tgt::clamp(current, 0, count()-1);

        if(event->modifiers() & Qt::Modifier::CTRL) {
            property_->moveInstance(instance.instanceId_, current);
        }
        else {
            property_->setSelectedInstance(current);
        }
    }

    void setPartner(QWidget* partner) {
        partner_ = partner;
    }

private:

    InteractiveListProperty* property_;
    QWidget* partner_;
    QPoint dragStartPosition_;
};


InteractiveListPropertyWidget::InteractiveListPropertyWidget(InteractiveListProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent)
    , itemTable_(nullptr)
    , instanceTable_(nullptr)
    , property_(prop)
{
    tgtAssert(property_, "null pointer");

    itemTable_ = new ItemListWidget(property_, nullptr);
    instanceTable_ = new InstanceListWidget(property_, nullptr);

    itemTable_->setToolTip(tr("Use Drag and Drop or press RETURN to add an instance of the selected item"));
    instanceTable_->setToolTip(tr("Use Drag and Drop or CTRL + UP/Down to move items"));

    // Reference partners.
    itemTable_->setPartner(instanceTable_);
    instanceTable_->setPartner(itemTable_);

    // Add widgets.
    addWidget(itemTable_);
    addWidget(instanceTable_);

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
