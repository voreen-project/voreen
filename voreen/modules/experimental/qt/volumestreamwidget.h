/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VOLUMESTREAMWIDGET_H
#define VOLUMESTREAMWIDGET_H

#include <tgt/types.h>
#include <QMutex>
#include <QThread>
#include <QWidget>
#include "tgt/types.h"

class QPushButton;

namespace voreen {

class VolumeStreamProperty;
class VolumeStreamProcessor;
class CanvasRenderer;

class VolumeStreamThread: public QThread {
public:
    VolumeStreamThread(VolumeStreamProcessor* proc, QMutex& worklistMutex);
    void run();

    bool running_;
    int runtime_;
protected:
    VolumeStreamProcessor* proc_;
    QMutex& worklistMutex_;
};

class VolumeStreamLoadThread: public VolumeStreamThread {
    Q_OBJECT
public:
    VolumeStreamLoadThread(VolumeStreamProcessor* proc, QMutex& worklistMutex);
    void run();
public slots:
    void loadStep(int step);
};

class VolumeStreamProcessThread: public VolumeStreamThread {
    Q_OBJECT
public:
    VolumeStreamProcessThread(VolumeStreamProcessor* proc, QMutex& worklistMutex);
    void run();
public slots:
    void processStep(int step);

signals:
    void stepProcessed(int step);
};

class VolumeStreamWidget : public QWidget {
    Q_OBJECT
public:
    VolumeStreamWidget(VolumeStreamProperty* prop, QWidget* parent = 0);
    ~VolumeStreamWidget();

public slots:
    void play();
    void showStep(int step);
    void stopThreads();

private:
    VolumeStreamProperty* prop_;
    VolumeStreamProcessor* proc_;

    QPushButton* playButton_;
    QMutex worklistMutex_;

    VolumeStreamLoadThread* loadThread_;
    VolumeStreamProcessThread* processThread_;
    VolumeStreamProcessThread* processThread2_;

    int timeInclude_;
    int timeRender_;

    CanvasRenderer* renderer_;

    uint64_t startTime_;
    int frameCounter_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMESTREAMWIDGET_H
