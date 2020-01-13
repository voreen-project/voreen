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

#include "bodyparts3dwidget.h"
#include "../datastructures/bodyparts3dmodel.h"

#include "tgt/filesystem.h"
#include "voreen/qt/voreenapplicationqt.h"

#include <QDialog>
#include <QColorDialog>
#include <QPushButton>
#include <QFileDialog>
#include <QGridLayout>
#include <QMainWindow>
#include <QMdiSubWindow>

namespace voreen {
using tgt::vec4;

const std::string BodyParts3DWidget::loggerCat_("voreen.qt.BodyParts3DWidget");

BodyParts3DWidget::BodyParts3DWidget(QWidget* parent, BodyParts3DSource* bodyParts3DSource)
    : QProcessorWidget(bodyParts3DSource, parent)
{
    mainLayout_ = new QVBoxLayout();

    model_ = new BodyParts3DModel(this);

    view_ = new QTreeView(this);
    view_->setModel(model_);
    view_->setContextMenuPolicy(Qt::CustomContextMenu);
    view_->setSelectionMode(QAbstractItemView::ExtendedSelection);


    saveAction_ = new QAction(QIcon(":/qt/icons/document_save.png"),"Save", this);
    loadAction_ = new QAction(QIcon(":/qt/icons/document_open.png"),"Load",this);
    loadBodyPartsDataAction_ = new QAction(QIcon(":/qt/icons/apply.png"),"Render", this);


    toolBar_ = new QToolBar("ToolBar", this);
    toolBar_->addAction(saveAction_);
    toolBar_->addAction(loadAction_);
    toolBar_->addSeparator();
    toolBar_->addAction(loadBodyPartsDataAction_);

    combo_ = new QComboBox(this);
    combo_->addItem("Version 3.0 - Conventional Part-of Relation");
    combo_->addItem("Version 3.0HQ - Conventional Part-of Relation");
    combo_->addItem("Version 4.0 - Is-A Tree");

    mainLayout_->addWidget(combo_);
    mainLayout_->addWidget(toolBar_);
    mainLayout_->addWidget(view_);

    setMinimumSize(711,400);
    setLayout(mainLayout_);

    tgtAssert(bodyParts3DSource, "No BodyParts3DSource");
    connect(combo_ , SIGNAL(currentIndexChanged(int)),this,SLOT(switchversion()));
    connect(saveAction_,SIGNAL(triggered()),this,SLOT(save()));
    connect(loadAction_,SIGNAL(triggered()),this,SLOT(load()));
    connect(loadBodyPartsDataAction_,SIGNAL(triggered()),this,SLOT(loadBodyPartsData()));
    connect(view_, SIGNAL(clicked(const QModelIndex&)), this, SLOT(clicked(const QModelIndex&)));
    path_ = "";
    version_ = Version3_0;
}

BodyParts3DWidget::~BodyParts3DWidget(){
    if(view_)
        delete view_;
    if(model_)
        delete model_;
    if(combo_)
        delete combo_;
}

void BodyParts3DWidget::updateFromProcessor(){
    BodyParts3DSource* proc = static_cast<BodyParts3DSource*>(processor_);
    view_->setModel(new BodyParts3DModel(this));
    delete model_;
    model_ = 0;
    if(proc){
        //if(proc->getPath() != path_ || proc->getVersion() != version_){
            proc->setPathValid(false);
            version_ = proc->getVersion();
            path_ = proc->getPath();
            if(version_ == Version3_0){
                combo_->setCurrentIndex(0);
                model_ = new BodyParts3DModel(path_, false, version_, this);
            }else if(version_ == Version3_0HQ){
                combo_->setCurrentIndex(1);
                model_ = new BodyParts3DModel(path_, false, version_, this);
            }else{
                combo_->setCurrentIndex(2);
                model_ = new BodyParts3DModel(path_, false, version_, this);
            }
            proc->setPathValid(true);
        //}
    }
    view_->setModel(model_);
}

void BodyParts3DWidget::switchversion(){
    if(combo_->currentText() == QString("Version 3.0 - Conventional Part-of Relation")){
        version_ = Version3_0;
    }else if(combo_->currentText() == QString("Version 3.0HQ - Conventional Part-of Relation")){
        version_ = Version3_0HQ;
    }else{
        version_ = Version4_0;
    }
}

void BodyParts3DWidget::save(){
    QString fileName = QFileDialog::getSaveFileName(this,tr("Save BodyParts3D ColorSheme"), "",tr("all files (*)"));
    model_->exportCatalog(fileName.toStdString());
}

void BodyParts3DWidget::load(){
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open BodyParts3D ColorSheme"), "", tr("all files (*)"));
    if(fileName != QString("")){
        path_ = fileName.toStdString();
        delete model_;
        model_ = 0;
        model_ = new BodyParts3DModel(path_, true, version_, this);
        path_ = path_.substr(0,path_.find_last_of('/'));
        view_->setModel(model_);
        view_->show();
        BodyParts3DSource* proc = static_cast<BodyParts3DSource*>(processor_);
        if(proc){
            proc->setVersion(version_);
            proc->setPath(path_);
        }
    }
}

void BodyParts3DWidget::clicked(const QModelIndex &index){
    if(index.column() == 2) {
        BodyParts3DItem* item = static_cast<BodyParts3DItem*>(index.internalPointer());
        if(item) {
            tgt::vec4 prevC = item->getColor();
            prevC *= 255.0f;
            QColor prevColor(prevC.x, prevC.y, prevC.z, prevC.w);
            QColor color = QColorDialog::getColor(prevColor, this);

            if (color.isValid()) {
                tgt::vec4 c(color.red(), color.green(), color.blue(), color.alpha());
                c /= 255.0f;
                item->setColor(c);
                processor_->invalidate();
            }
        }
    }
}

void BodyParts3DWidget::loadBodyPartsData(){
    std::vector<std::string> fileNameList = model_->getFileNames();
    std::vector<tgt::vec4> colorList = model_->getColors();
    std::vector<std::string> descList = model_->getDescriptions();
    BodyParts3DSource *proc = static_cast<BodyParts3DSource*>(processor_);
    if(proc){
        proc->createMesh(fileNameList, colorList, descList);
    }
}


}// namespace
