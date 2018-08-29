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

#ifndef VRN_VOREENQTSTARTUPWIZARD_H
#define VRN_VOREENQTSTARTUPWIZARD_H
#include "voreen/qt/voreenqtapi.h"

#include <QDialog>
#include <QLabel>
#include <QGridLayout>
#include <QTabWidget>
#include <QListWidget>
#include <QPushButton>

namespace voreen{
class VRN_QT_API VoreenStartupWizard :public QDialog{
        Q_OBJECT
public:
    VoreenStartupWizard( QStringList recentWorkspaceNames, QStringList standardWorkspaceNames, QWidget *parent = 0, Qt::WindowFlags f = 0);
    ~VoreenStartupWizard();
    QString getSelectedWorkspace();
    virtual QLabel* getLogoLabel();

protected:
    QStringList recentWorkspaceNames_;
    QStringList standardWorkspaceNames_;

    QBoxLayout* windowsLayout_;

    QTabWidget* tabWidget_;
    QListWidget* recentWorkspaces_;
    QListWidget* stadardWorkspaces_;

    QBoxLayout* buttonLayout_;
    QPushButton* openButton_;
    QPushButton* openWorkspaceFromFileButton_;
    //QPushButton* newWorkspaceButton_;

    QLabel* logo_;

    QString forceSelected_;

    QListWidget* createWorkspaceListWidget(QStringList list);
private slots:
    void openItem(QListWidgetItem*);
    void currentTabChanged(int index);
    void openWorkspaceFromFile(void);
};
};

#endif
