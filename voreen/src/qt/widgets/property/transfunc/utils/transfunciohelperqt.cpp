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

#include "voreen/qt/widgets/property/transfunc/utils/transfunciohelperqt.h"

#include "voreen/core/voreenapplication.h"
//tf
#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include "tgt/logmanager.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QList>
#include <QUrl>
#include <QDesktopServices>

namespace voreen {

const QString TransFuncIOHelperQt::getLoadFileName(QStringList filters) {
    //create dialog
    QFileDialog fileDialog(0,tr("Choose a transfer function to open"),VoreenApplication::app()->getCoreResourcePath("transferfuncs").c_str());
    fileDialog.setNameFilters(filters);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getCoreResourcePath("transferfuncs").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    if (fileDialog.exec() && !fileDialog.selectedFiles().empty()) {
        return fileDialog.selectedFiles()[0];
    }

    return QString();
}

const QString TransFuncIOHelperQt::getSaveFileName(QStringList filters) {
    //create dialog
    QFileDialog fileDialog(0, tr("Choose a filename to save transfer function"), VoreenApplication::app()->getCoreResourcePath("transferfuncs").c_str());
    fileDialog.setNameFilters(filters);
    fileDialog.setAcceptMode(QFileDialog::AcceptSave);
    fileDialog.setOption(QFileDialog::DontUseNativeDialog);

    QList<QUrl> urls;
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getCoreResourcePath("transferfuncs").c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getUserDataPath().c_str());
    urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("modules").c_str());
    if (QDir(VoreenApplication::app()->getBasePath("custommodules").c_str()).exists())
        urls << QUrl::fromLocalFile(VoreenApplication::app()->getBasePath("custommodules").c_str());
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::DesktopLocation))
        urls << QUrl::fromLocalFile(f);
    for (auto f : QStandardPaths::standardLocations(QStandardPaths::HomeLocation))
        urls << QUrl::fromLocalFile(f);
    fileDialog.setSidebarUrls(urls);

    QStringList fileList;
    if (fileDialog.exec() && !fileDialog.selectedFiles().empty()) {
        QString endingFilter = fileDialog.selectedNameFilter();
        int pos = endingFilter.lastIndexOf(".");
        //removes closing bracket
        endingFilter.chop(1);
        endingFilter = endingFilter.mid(pos);

        //look whether the user specified an ending
        fileList = fileDialog.selectedFiles();
        size_t dotPosition = fileList[0].toStdString().rfind(".");
        if (dotPosition == std::string::npos) {
            // no ending given -> add ending of selected filter
            fileList[0].append(endingFilter);
        }
        else {
            // an ending was given -> test whether it matches the selected filter
            if (fileList[0].mid(static_cast<int>(dotPosition)) != endingFilter)
                fileList[0].append(endingFilter);
        }
        return fileList[0];
    }
    return QString();
}


bool TransFuncIOHelperQt::loadTransferFunction(TransFuncBase* func) {
    tgtAssert(func,"empty function");
    //create filter with supported file formats
    QStringList filters;
    std::string temp = "transfer function (";
    for (size_t i = 0; i < func->getLoadFileFormats().size(); ++i) {
        temp += "*." + func->getLoadFileFormats()[i] + " ";
    }
    temp.replace(temp.length()-1,1,")");
    filters << temp.c_str();
    for (size_t i = 0; i < func->getLoadFileFormats().size(); ++i) {
        std::string temp = "transfer function (*." + func->getLoadFileFormats()[i] + ")";
        filters << temp.c_str();
    }

    QString fileName = getLoadFileName(filters);
    if (!fileName.isEmpty()) {
        if (!func->load(fileName.toStdString())) {
            QMessageBox::critical(0, tr("Error"),
                "The selected transfer function could not be loaded.");
            LERRORC("TrunsFuncIOHelperQt","The selected transfer function could not be loaded. Maybe the file is corrupt.");
            return false;
        }
        return true;
    }
    return false;
}

bool TransFuncIOHelperQt::saveTransferFunction(TransFuncBase* func) {
    tgtAssert(func,"empty function");
    //create filter with supported file formats
    QStringList filters;
    for (size_t i = 0; i < func->getSaveFileFormats().size(); ++i) {
        std::string temp = "transfer function (*." + func->getSaveFileFormats()[i] + ")";
        filters << temp.c_str();
    }

    QString fileName = getSaveFileName(filters);
    if (!fileName.isEmpty()) {
        //save transfer function to disk
        if (!func->save(fileName.toStdString())) {
            QMessageBox::critical(0, tr("Error"),
                                  tr("The transfer function could not be saved."));
            LERRORC("TransFuncIOHelperQt","The transfer function could not be saved. Maybe the disk is full?");
            return false;
        }
        return true;
    }
    return false;
}



} // namespace voreen
