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

#ifndef VRN_PROGRESSDIALOG_H
#define VRN_PROGRESSDIALOG_H

#include <QProgressDialog>
#include <QTime>
#include <memory>

#include "voreen/core/io/progressbar.h"

namespace voreen {

/**
 * A class for a small dialog indicating the progress of loading a file by
 * using a QProgressDialog.
 * Note: title and message, both inherited by ProgressBar, need to be set before
 *       the first update is triggered.
 *       Update will be triggered by a successful call of setProgress() or show().
 */
class ProgressDialog : public ProgressBar {
public:
    ProgressDialog(QWidget* parent);
    ~ProgressDialog();

    void update();
    void show();
    void hide();
    void forceUpdate();

protected:

    QWidget* parent_;
    QTime updateTime_;
    std::unique_ptr<QProgressDialog> progressDialog_;
};

} // namespace

#endif // VRN_PROGRESSDIALOG_H
