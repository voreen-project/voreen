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

#include "labelingwidgetqt.h"

#include "voreen/core/utils/voreenpainter.h"

#include <QWidget>
#include <QFrame>
#include <QPushButton>

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QInputDialog>

namespace voreen {

LabelingWidgetQt::LabelingWidgetQt(QWidget *parent) :
    QWidget(parent)
{
    setObjectName(tr("Labeling"));
    createWidgets();
    createConnections();
}

void LabelingWidgetQt::createWidgets() {
    // resize(300,300);
    QVBoxLayout* mainLayout = new QVBoxLayout;

    tabWidget_ = new QTabWidget();
    mainLayout->addWidget(tabWidget_);

    QListWidget* heartCTList = new QListWidget();
    heartCTList->setSelectionMode(QAbstractItemView::MultiSelection);
    heartCTList->addItem("Columna vertebralis");
    heartCTList->addItem("Aorta");
    heartCTList->addItem("V. cava inferior");
    heartCTList->addItem("V. cava superior");
    heartCTList->addItem("Atrium dextrum");
    heartCTList->addItem("Ventriculus sinister");
    heartCTList->addItem("Ventriculus dexter");
    heartCTList->addItem("Truncus pulmonalis");
    heartCTList->addItem("A. pulmonalis sinistra");
    heartCTList->addItem("A. pulmonalis dextra");
    heartCTList->addItem("Atrium sinistrum");
    heartCTList->addItem("V. pulmonales sinistrae");
    heartCTList->addItem("V. pulmonales dextrae");
    tabWidget_->addTab(heartCTList, "Heart CT");
    listWidgets_.push_back(heartCTList);

    QListWidget* bodyCTList = new QListWidget();
    bodyCTList->setSelectionMode(QAbstractItemView::MultiSelection);
    bodyCTList->addItem("Item1");
    bodyCTList->addItem("Item2");
    tabWidget_->addTab(bodyCTList, "Body CT");
    listWidgets_.push_back(bodyCTList);

    QListWidget* skullMRTList = new QListWidget();
    skullMRTList->setSelectionMode(QAbstractItemView::MultiSelection);
    skullMRTList->addItem("Item1");
    skullMRTList->addItem("Item2");
    tabWidget_->addTab(skullMRTList, "Skull MRT");
    listWidgets_.push_back(skullMRTList);

    QListWidget* handCTList = new QListWidget();
    handCTList->setSelectionMode(QAbstractItemView::MultiSelection);
    handCTList->addItem("Item1");
    handCTList->addItem("Item2");
    tabWidget_->addTab(handCTList, "Hand CT");
    listWidgets_.push_back(handCTList);

    QHBoxLayout* hboxLayout1 = new QHBoxLayout();
    hboxLayout1->addWidget(new QPushButton("Add All"));
    hboxLayout1->addWidget(new QPushButton("Remove All"));
    mainLayout->addItem(hboxLayout1);
    QFrame* separator = new QFrame();
    separator->setFrameShape(QFrame::HLine);
    separator->setFrameShadow(QFrame::Sunken);
    mainLayout->addWidget(separator);
    QHBoxLayout* hboxLayout2 = new QHBoxLayout();
    hboxLayout2->addWidget(new QPushButton("New Label"));
    hboxLayout2->addWidget(new QPushButton("Edit Label"));
    hboxLayout2->addWidget(new QPushButton("Delete Label"));
    mainLayout->addItem(hboxLayout2);

    setLayout(mainLayout);
}

void LabelingWidgetQt::createConnections() {
    for (size_t i=0; i<listWidgets_.size(); ++i){
        //connect( listWidgets_[i], SIGNAL(itemSelectionChanged()), this, SLOT(selectionChanged()) );
        connect( listWidgets_[i], SIGNAL(itemPressed(QListWidgetItem*)), this, SLOT(selectionChanged(QListWidgetItem*)) );
    }
}

void LabelingWidgetQt::selectionChanged(QListWidgetItem* /*item*/) {
    //std::cout << "selection changed" << std::endl;
    /*if (item->listWidget()->selectedItems().contains(item))
        postMessage(new StringMsg("add.unsegmentedLabel", item->text().toStdString() ) );
    else
        postMessage(new StringMsg("remove.unsegmentedLabel", item->text().toStdString() ) );
    */
}

std::string LabelingWidgetQt::showInputDialog(const std::string &title,
                                                   const std::string &label,
                                                   const std::string &text)
{
    return QInputDialog::getText(NULL, QString::fromStdString(title),
                QString::fromStdString(label),
                QLineEdit::Normal,
                QString::fromStdString(text) ).toStdString();
}

}  // namespace
