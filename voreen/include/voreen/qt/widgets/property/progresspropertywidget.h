/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_PROGRESSPROPERTYWIDGET_H
#define VRN_PROGRESSPROPERTYWIDGET_H

#include "voreen/qt/widgets/property/qpropertywidget.h"

#include <QTime>

class QProgressBar;

namespace voreen {

class ProgressProperty;

class ProgressPropertyWidget : public QPropertyWidget {
Q_OBJECT
public:
    ProgressPropertyWidget(ProgressProperty* prop, QWidget* parent = 0);

protected slots:
    virtual void updateFromPropertySlot();

private:
    QTime time_;    ///< prevent crash in release mode
    ProgressProperty* property_;
    QProgressBar* progressBar_;

};

} // namespace

#endif // VRN_BUTTONPROPERTYWIDGET_H
