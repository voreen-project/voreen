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

#include "voreenbiologystartupwizard.h"
namespace voreen{
VoreenBiologyStartupWizard::VoreenBiologyStartupWizard(QStringList recentWorkspaceNames, QStringList standardWorkspaceNames, 
                                             QWidget *parent, Qt::WindowFlags f)
    :VoreenStartupWizard(recentWorkspaceNames, standardWorkspaceNames, parent, f){
    // We need to do this here, because in the constructor of
    // VoreenStartupWizard the correct v-table-pointer is not set
    // yet.
    windowsLayout_->addWidget(getLogoLabel());
}
    
QLabel* VoreenBiologyStartupWizard::getLogoLabel(){
    QPixmap logoPixmap(":/voreenbiology/image/logo-startup-wizard.png");
    QLabel *logo = new QLabel();
    logo->setPixmap(logoPixmap);
    return logo;
}
}