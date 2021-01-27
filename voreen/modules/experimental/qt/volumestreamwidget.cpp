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

#include "volumestreamwidget.h"

#include "../processors/volumestreamprocessor.h"

#include "voreen/core/voreenapplication.h"
#include "../properties/volumestreamproperty.h"

#include "voreen/core/network/processornetwork.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/processors/processorwidget.h"

// core module is always available
#include "modules/core/processors/output/canvasrenderer.h"
#include "modules/core/qt/processor/canvasrendererwidget.h"

#include <QApplication>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>

namespace voreen {

typedef VolumeStreamProcessor::StreamStep StreamStep;

const std::string VolumeStreamWidget::loggerCat_("voreen.VolumeStreamWidget");

namespace {
int const MAX_UNFINISHED_STEPS = 5;

uint64_t getTicks() {
    uint64_t ticks;

#ifdef __unix__
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    ticks = now.tv_sec * 1000000 + now.tv_nsec / 1000;
#else
    ticks = 0;
#endif
    return ticks;
}

} // namespace

VolumeStreamThread::VolumeStreamThread(VolumeStreamProcessor* proc, QMutex& worklistMutex)
    : running_(false),
      runtime_(0),
      proc_(proc),
      worklistMutex_(worklistMutex)
{
}

void VolumeStreamThread::run() {
    exec();
}

VolumeStreamLoadThread::VolumeStreamLoadThread(VolumeStreamProcessor* proc, QMutex& worklistMutex)
    : VolumeStreamThread(proc, worklistMutex)
{}

void VolumeStreamLoadThread::run() {
    std::deque<StreamStep>& worklist = proc_->getWorkList();
    std::deque<StreamStep>::iterator it;
    running_ = true;
    while (running_) {
        int step = -1;
        worklistMutex_.lock();
        if (worklist.empty()) {
            worklistMutex_.unlock();
            return;
        }

        it = worklist.begin();
        int unfinished = 0;
        while (it != worklist.end()) {
            if (!(*it).loaded_) {
                step = (*it).step_;
                break;
            }
            ++it;
            unfinished++;
        }

        worklistMutex_.unlock();

        if (step >= 0 && unfinished < MAX_UNFINISHED_STEPS)
            loadStep(step);
        else
            msleep(2);
    }
}

void VolumeStreamLoadThread::loadStep(int step) {
    std::deque<StreamStep>& worklist = proc_->getWorkList();
    std::deque<StreamStep>::iterator it;

    worklistMutex_.lock();

    StreamStep item;
    it = worklist.begin();
    while (it != worklist.end()) {
        if ((*it).step_ == step) {
            item = *it;
            break;
        }
        ++it;
    }

    if (item.loaded_) {
        worklistMutex_.unlock();
        return;
    }

    worklistMutex_.unlock();

    // load
    uint64_t tick = getTicks();
    proc_->readBlocks(step);
    int time = (getTicks() - tick) / 1000;
    runtime_ += time;
    // std::cout << "load " << step << " (" << time << ")" << std::endl;
    worklistMutex_.lock();

    it = worklist.begin();
    while (it != worklist.end()) {
        if ((*it).step_ == item.step_ && !(*it).loaded_) {
            (*it).data_ = item.data_;
            (*it).loaded_ = true;
            break;
        }
        ++it;
    }

    worklistMutex_.unlock();
}

VolumeStreamProcessThread::VolumeStreamProcessThread(VolumeStreamProcessor* proc, QMutex& worklistMutex)
    : VolumeStreamThread(proc, worklistMutex)
{}

void VolumeStreamProcessThread::run() {
    std::deque<StreamStep>& worklist = proc_->getWorkList();
    std::deque<StreamStep>::iterator it;
    running_ = true;
    while (running_) {
        int step = -1;
        worklistMutex_.lock();
        if (worklist.empty()) {
            worklistMutex_.unlock();
            return;
        }

        it = worklist.begin();
        while (it != worklist.end()) {
            if ((*it).loaded_ && !(*it).processed_) {
                step = (*it).step_;
                break;
            }
            ++it;
        }

        worklistMutex_.unlock();

        if (step >= 0)
            processStep(step);
        else
            msleep(2);
    }
}

void VolumeStreamProcessThread::processStep(int step) {
    std::deque<StreamStep>& worklist = proc_->getWorkList();
    std::deque<StreamStep>::iterator it;

    worklistMutex_.lock();

    StreamStep item;
    it = worklist.begin();
    while (it != worklist.end()) {
        if ((*it).step_ == step) {
            item = *it;
            break;
        }
        ++it;
    }

    if (item.processed_) {
        worklistMutex_.unlock();
        return;
    }

    worklistMutex_.unlock();

    // process
    uint64_t tick = getTicks();
    proc_->uncompressBlocks(step);
    int time = (getTicks() - tick) / 1000;
    runtime_ += time;
    // std::cout << "proc " << step << " (" << time << ")" << std::endl;

    worklistMutex_.lock();

    it = worklist.begin();
    while (it != worklist.end()) {
        if ((*it).step_ == item.step_ && !(*it).processed_) {
            (*it).data_ = item.data_;
            (*it).processed_ = true;
            break;
        }
        ++it;
    }

    worklistMutex_.unlock();

    emit stepProcessed(step);
}

VolumeStreamWidget::VolumeStreamWidget(VolumeStreamProperty* prop, QWidget* parent)
    : QWidget(parent),
      prop_(prop),
      proc_(dynamic_cast<VolumeStreamProcessor*>(prop->getOwner())),
      timeInclude_(0),
      timeRender_(0),
      renderer_(0),
      startTime_(0)
{
    setObjectName("VolumeStreamWidget");

    if (!proc_)
        LERROR("Could not get VolumeStreamProcessor via property owner!");

    loadThread_ = 0;
    processThread_ = 0;
    processThread2_ = 0;


    QVBoxLayout* mainLayout = new QVBoxLayout();

    // QGroupBox* orientationBox_ = new QGroupBox(tr("Orientation and Distance"));
    // QGridLayout* gridLayout = new QGridLayout();
    // gridLayout->setColumnStretch(0, 3);
    // gridLayout->setColumnStretch(1, 2);
    // gridLayout->addWidget(new QLabel(tr("Orientation: ")), 1, 0);
    // orientationBox_->setLayout(gridLayout);
    // mainLayout->addWidget(orientationBox_);


    playButton_ = new QPushButton(tr("Play"));
    mainLayout->addWidget(playButton_);

    mainLayout->addStretch();


    setLayout(mainLayout);

    connect(playButton_, SIGNAL(clicked()), this, SLOT(play()));
}

VolumeStreamWidget::~VolumeStreamWidget() {
}

void VolumeStreamWidget::play() {
    const std::vector<Processor*> processors
        = VoreenApplication::app()->getNetworkEvaluator()->getProcessorNetwork()->getProcessors();
    for (std::vector<Processor*>::const_iterator iter = processors.begin(); iter != processors.end(); iter++) {
        CanvasRenderer* canvasProc = dynamic_cast<CanvasRenderer*>(*iter);
        if (canvasProc) {
            renderer_ = canvasProc;
            break;
        }
    }

    if (loadThread_) {
        stopThreads();
        return;
    }

    proc_->flush();
    proc_->setStep(0);
    proc_->loadStep();

    frameCounter_ = 0;

    if (proc_->getBlockSize() == 0) {
        startTime_ = getTicks();
        for (int i=0; i <= proc_->getStepCount(); i++) {
            proc_->setStep(i);
            if (renderer_)
                renderer_->getCanvas()->repaint();
            frameCounter_++;
        }
        uint64_t time = getTicks() - startTime_;
        std::cout << "TOTAL time: " << (time / 1000000.f) << " sec, "
                  << frameCounter_ / (time / 1000000.f) << " fps" << std::endl;
        std::cout << "TOTAL disk: " << proc_->getRawSize() / (1024*1024) << " MB, "
                  << (proc_->getRawSize() / (1024.f*1024.f)) / (time / 1000000.f) << " MB/s" << std::endl;
        return;
    }

    loadThread_ = new VolumeStreamLoadThread(proc_, worklistMutex_);
    processThread_ = new VolumeStreamProcessThread(proc_, worklistMutex_);
//    processThread2_ = new VolumeStreamProcessThread(proc_, worklistMutex_);


    std::deque<StreamStep>& worklist = proc_->getWorkList();
    worklist.clear();

    connect(processThread_, SIGNAL(stepProcessed(int)),
            this, SLOT(showStep(int)), Qt::QueuedConnection);
    // connect(processThread2_, SIGNAL(stepProcessed(int)),
    //         this, SLOT(showStep(int)), Qt::QueuedConnection);

    for (int i=0; i <= proc_->getStepCount(); i++) {
        worklistMutex_.lock();
        worklist.push_back(i);
        worklistMutex_.unlock();
    }

    startTime_ = getTicks();
    timeInclude_ = 0;
    timeRender_ = 0;

    loadThread_->start();
    processThread_->start();
//    processThread2_->start();
}

void VolumeStreamWidget::showStep(int step) {
    const std::vector<Processor*> processors
        = VoreenApplication::app()->getNetworkEvaluator()->getProcessorNetwork()->getProcessors();
    for (std::vector<Processor*>::const_iterator iter = processors.begin(); iter != processors.end(); iter++) {
        CanvasRenderer* canvasProc = dynamic_cast<CanvasRenderer*>(*iter);
        if (canvasProc) {
            renderer_ = canvasProc;
            break;
        }
    }

    worklistMutex_.lock();

    std::deque<StreamStep>& worklist = proc_->getWorkList();
    StreamStep item;
    std::deque<StreamStep>::iterator it = worklist.begin();
    while (it != worklist.end()) {
        if ((*it).step_ == step) {
            worklist.erase(it);
            break;
        }
        ++it;
    }

    worklistMutex_.unlock();

    uint64_t tick = getTicks();

    // TODO: valid code?
    if (step > 0)
        proc_->includeBlocks(step, 1);
    int timeIncl = (getTicks() - tick) / 100;
    timeInclude_ += timeIncl;
//    std::cout << "incl " << step << " (" << timeIncl << ")" << std::endl;

    tick = getTicks();

    proc_->showStep(step);
    // We must call repaint() directly here, instead of going through the Qt event queue,
    // because this method is called as a slot from a different thread.
    if (renderer_)
        renderer_->getCanvas()->repaint();
    int timeRend = (getTicks() - tick) / 1000;
    timeRender_ += timeRend;

//    std::cout << "rend " << step << " (" << timeRend << ")" << std::endl;

    frameCounter_++;

    if (worklist.empty())
        stopThreads();
}

void VolumeStreamWidget::stopThreads() {
    if (loadThread_) {
        loadThread_->running_ = false;
        loadThread_->wait();

        std::cout << "TOTAL load: " << loadThread_->runtime_ << std::endl;
        delete loadThread_;
        loadThread_ = 0;
    }

    if (processThread_) {
        processThread_->running_ = false;
        processThread_->wait();

        std::cout << "TOTAL proc: " << processThread_->runtime_ << std::endl;
        delete processThread_;
        processThread_ = 0;
    }

    if (processThread2_) {
        processThread2_->running_ = false;
        processThread2_->wait();

        std::cout << "TOTAL proc2: " << processThread2_->runtime_ << std::endl;
        delete processThread2_;
        processThread2_ = 0;
    }

    std::cout << "TOTAL incl: " << timeInclude_ << std::endl;
    std::cout << "TOTAL rend: " << timeRender_ << std::endl;

    uint64_t time = getTicks() - startTime_;

    std::cout << "TOTAL time: " << (time / 1000000.f) << " sec, "
              << frameCounter_ / (time / 1000000.f) << " fps" << std::endl;
    std::cout << "TOTAL disk: " << proc_->getRawSize() / (1024*1024) << " MB, "
              << (proc_->getRawSize() / (1024.f*1024.f)) / (time / 1000000.f) << " MB/s" << std::endl;
    std::cout << "TOTAL step: " << frameCounter_ << std::endl;
}

} // namespace voreen
