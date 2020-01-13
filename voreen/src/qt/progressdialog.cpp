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

#include "voreen/qt/progressdialog.h"

#include <tgt/logmanager.h>
#include <tgt/assert.h>

#include <QCoreApplication>
#include <QThread>

namespace {
    const int stepGranularity = 200; // number of steps on the QProgressDialog
    const int minimalUpdateWait = 50; // time (ms) to wait between progress bar updates
    const int minimalShowWait = 400; // time (ms) to wait until the dialog appears
}

namespace voreen {

ProgressDialog::ProgressDialog(QWidget* parent)
    : parent_(parent)
{
}

ProgressDialog::~ProgressDialog() {
}

void ProgressDialog::update() {
    // GUI operations are only allowed in the GUI thread
    if (QThread::currentThread() != QCoreApplication::instance()->thread()) {
        tgtAssert(false, "GUI thread was not active, update discarded");
        LWARNINGC("qt.ProgressDialog", "GUI thread was not active, update discarded");
        return;
    }

    int intProgress = static_cast<int>(progress_ * stepGranularity);

    // Create the dialog as soon as the first update is done.
    // This is either done by show() or simply by setProgress();
    // It will also be freshly created if hide() was called before.
    if (!progressDialog_) {
        if (updateTime_.isNull()) {
            updateTime_.start();
        }
        else if(updateTime_.elapsed() >= minimalShowWait) {
            progressDialog_.reset(new QProgressDialog(parent_));
            progressDialog_->setCancelButton(nullptr);
            progressDialog_->setWindowTitle(QString::fromStdString(getTitle()));
            progressDialog_->setLabelText(QString::fromStdString(getMessage()));
            progressDialog_->setWindowModality(Qt::WindowModal);
            progressDialog_->setRange(0, stepGranularity);
            progressDialog_->setMinimumDuration(0);
            progressDialog_->setMinimumWidth(300);
            progressDialog_->setAutoReset(false);
            progressDialog_->setAutoClose(false);
            progressDialog_->setValue(intProgress);
            progressDialog_->show();
            progressDialog_->raise();
            progressDialog_->activateWindow();
        }
    }
    else if ((updateTime_.elapsed() >= minimalUpdateWait || intProgress == progressDialog_->maximum())
        && progressDialog_->value() != intProgress)
    {
        progressDialog_->setValue(intProgress);
        updateTime_.restart();
    }
}

void ProgressDialog::show() {
    // Trigger an update.
    setProgress(0.f);
}

void ProgressDialog::hide() {
    // Reset progress.
    setProgress(1.0f);
    // Delete the dialog.
    progressDialog_.reset();
    // Reset time to null.
    updateTime_ = QTime();
}

void ProgressDialog::forceUpdate() {
    // GUI operations are only allowed in the GUI thread
    if (QThread::currentThread() != QCoreApplication::instance()->thread()) {
        tgtAssert(false, "GUI thread was not active, update discarded");
        LWARNINGC("qt.ProgressDialog", "GUI thread was not active, update discarded");
        return;
    }
    
    if(!progressDialog_)
        update();
    else
        progressDialog_->repaint();
}

} // namespace
