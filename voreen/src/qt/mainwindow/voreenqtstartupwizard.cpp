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

#include "voreen/qt/mainwindow/voreenqtstartupwizard.h"
#include "voreen/core/voreenapplication.h"
#include <QFileInfo>
#include <QApplication>
#include <QFileDialog>
#include <QUrl>
#include <QDesktopServices>
#include "tgt/assert.h"
namespace voreen{

namespace{
class ListWidgetItem : public QListWidgetItem
{
public:
    ListWidgetItem(QString path);
    QString getPath();
private:
    QString path_;
};

ListWidgetItem::ListWidgetItem(QString path)
    :path_(path)
{
        QFileInfo fileInfo(path);
        setText(fileInfo.fileName());
        setToolTip(fileInfo.absoluteFilePath());
        if (!fileInfo.exists()){
            setForeground(Qt::red);
        }
}

QString ListWidgetItem::getPath(){
    return path_;
}

}


VoreenStartupWizard::VoreenStartupWizard( QStringList recentWorkspaceNames, QStringList stadardWorkspaceNames, QWidget *parent, Qt::WindowFlags f)
    : QDialog(parent, f)
    , recentWorkspaceNames_(recentWorkspaceNames)
    , standardWorkspaceNames_(stadardWorkspaceNames)
{

    buttonLayout_ = new QHBoxLayout();
    windowsLayout_ = new QVBoxLayout();

    tabWidget_ = new QTabWidget();

    recentWorkspaces_ = createWorkspaceListWidget(recentWorkspaceNames_);
    stadardWorkspaces_ = createWorkspaceListWidget(standardWorkspaceNames_);
    this->setContentsMargins(0, 0, 0, 0);
    windowsLayout_->setContentsMargins(0, 0, 0, 0);



    openButton_ = new QPushButton(tr("Load selected workspace"));
    openButton_->setDefault(true);
    connect(openButton_, SIGNAL(clicked()), this, SLOT(accept()));
    //newWorkspaceButton_ = new QPushButton(tr("New workspace"));
    //connect(newWorkspaceButton_, SIGNAL(clicked()), this, SLOT(reject()));
    openWorkspaceFromFileButton_  = new QPushButton(tr("Load Workspace from file"));
    connect(openWorkspaceFromFileButton_, SIGNAL(clicked()), this, SLOT(openWorkspaceFromFile()));



    tabWidget_->insertTab(0, recentWorkspaces_, "Recent workspaces");
    tabWidget_->insertTab(1, stadardWorkspaces_, "Template workspaces");
    connect(tabWidget_, SIGNAL(currentChanged(int)), this, SLOT(currentTabChanged(int)));

    if (recentWorkspaces_->count() != 0){
        tabWidget_->setCurrentIndex(0);
        // we are in tab 0 by default, so qt does
        // not call the slot that does initialization
        currentTabChanged(0);
    }else{
        // no recent workspaces exist
        // select template Workspace tab
        tabWidget_->setCurrentIndex(1);
    }


    setLayout(windowsLayout_);
    windowsLayout_->setSpacing(0);
    windowsLayout_->addWidget(tabWidget_);
    windowsLayout_->addSpacing(10);
    windowsLayout_->addLayout(buttonLayout_);
    windowsLayout_->addSpacing(10);

    buttonLayout_->setSpacing(10);
    buttonLayout_->addWidget(openWorkspaceFromFileButton_, 0, Qt::AlignLeft);
    buttonLayout_->addStretch();
    //buttonLayout_->addWidget(newWorkspaceButton_, 0, Qt::AlignRight);
    buttonLayout_->addWidget(openButton_, 0, Qt::AlignRight);
}
VoreenStartupWizard::~VoreenStartupWizard(){
    recentWorkspaces_->reset();
    stadardWorkspaces_->reset();
}

QString VoreenStartupWizard::getSelectedWorkspace(){
    if (!forceSelected_.isEmpty()){
        return forceSelected_;
    }

    QListWidget* widget = dynamic_cast<QListWidget*>(tabWidget_->currentWidget());

    ListWidgetItem *item = dynamic_cast<ListWidgetItem*>(widget->currentItem());
    return item->getPath();
}

QListWidget* VoreenStartupWizard::createWorkspaceListWidget(QStringList list){

    QListWidget* listWidget  = new QListWidget(this);
    for(int i = 0; i != list.size(); i++){
        QListWidgetItem* item = new ListWidgetItem(list[i]);

        listWidget->addItem(item);
    }

    listWidget->setCurrentRow(0);
    listWidget->setSelectionMode(QAbstractItemView::SingleSelection);
    listWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    listWidget->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    connect(listWidget, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(openItem(QListWidgetItem*)));
    return listWidget;
}

void VoreenStartupWizard::openItem(QListWidgetItem* qitem){
    ListWidgetItem *item = dynamic_cast<ListWidgetItem*>(qitem);
    forceSelected_ = item->getPath();
    emit accept();
}

QLabel* VoreenStartupWizard::getLogoLabel(){
    return new QLabel;
}

void VoreenStartupWizard::currentTabChanged(int index){
    QListWidget *listWidget = static_cast<QListWidget*>(tabWidget_->widget(index));

    listWidget->setFocus();

    // disable open button if there is nothing to be opend
    bool hasItems = listWidget->count() != 0;
    openButton_->setEnabled(hasItems);
}

void VoreenStartupWizard::openWorkspaceFromFile(void){
    QFileDialog fileDialog(this, tr("Open Workspace..."));

    QStringList filters;
    filters << "Voreen workspaces (*.vws)";
    fileDialog.setNameFilters(filters);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getApplicationResourcePath("workspaces").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    for (auto& f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto& f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec()){
        forceSelected_ = fileDialog.selectedFiles().at(0);
        emit accept();
    }
}

}
