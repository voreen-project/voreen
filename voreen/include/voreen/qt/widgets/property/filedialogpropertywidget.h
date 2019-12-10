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

#ifndef VRN_FILEDIALOGPROPERTYWIDGET_H
#define VRN_FILEDIALOGPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

class QPushButton;
class QCheckBox;
class QVBoxLayout;

namespace voreen {

class FileDialogProperty;

class FileDialogPropertyWidget : public QPropertyWidget {
    Q_OBJECT
public:
    FileDialogPropertyWidget(FileDialogProperty* prop, QWidget* parent = 0);
    virtual ~FileDialogPropertyWidget();

public slots:
    void setProperty();
    void changeFileWatchOption(bool checked);

protected:
    void updateButtonText(const std::string& filename);

    FileDialogProperty* property_; 
    QPushButton* openFileDialogBtn_;
    QCheckBox* fileWatchCheckBox_;
    QVBoxLayout* mainLayout_;

protected slots:
    virtual void updateFromPropertySlot();

};

} // namespace

#endif // VRN_FILEDIALOGPROPERTYWIDGET_H
