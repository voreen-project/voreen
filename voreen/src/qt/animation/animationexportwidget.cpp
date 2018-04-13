/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/qt/animation/animationexportwidget.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/networkevaluator.h"

#include "modules/core/processors/output/canvasrenderer.h" //< core module is always available

#include "tgt/qt/qtcanvas.h"

#include <QApplication>
#include <QLabel>
#include <QMessageBox>
#include <QProgressDialog>
#include <QTimeEdit>
#include <QVariant>

#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

AnimationExportWidget::AnimationExportWidget(QWidget* parent, Animation* animation, NetworkEvaluator* network,
                                             float duration, int startframe, int endframe, float stretchFactor)
    : QDialog(parent)
    , animation_(animation)
    , network_(network)
    , duration_(duration)
    , startframe_(startframe)
    , endframe_(endframe)
    , canvasRenderer_(0)
    , canvasSize_(256, 256)
    , renderingVideo_(false)
    , recordPathName_("")
    , fpsFactor_(stretchFactor)
    , renderState_(Inactive)
{
    setWindowTitle(tr("Export Animation"));
    setObjectName(tr("Export Options"));
    createWidgets();
    createConnections();
}

AnimationExportWidget::~AnimationExportWidget(){
    comboCanvases_->disconnect();
    renderBox_->disconnect();
}
void AnimationExportWidget::networkChanged() {
    if (network_ == 0)
        return;

    canvasRenderer_ = 0;
    allCanvasRenderer_.clear();

    const std::vector<Processor*>& processors = network_->getProcessorNetwork()->getProcessors();
    for (size_t i = 0; i < processors.size(); ++i) {
        if (CanvasRenderer* cr = dynamic_cast<CanvasRenderer*>(processors[i])) {
            allCanvasRenderer_.push_back(cr);
        }
    }

    if (!allCanvasRenderer_.empty()) {
        canvasRenderer_ = allCanvasRenderer_[0];
        canvasSize_ = canvasRenderer_->getCanvasSize();
        int width = canvasSize_.x;
        if (width % 2 == 1)
            width -= 1;

        int height = canvasSize_.y;
        if (height % 2 == 1)
            height -= 1;

        spinWidth_->setValue(width);
        spinHeight_->setValue(height);
    }

    refreshComboBoxes();
}

void AnimationExportWidget::controlledCanvasChanged(int index) {
    if ((comboCanvases_ == 0) || (index < 0) || (index >= comboCanvases_->count()))
        return;
    canvasRenderer_ = static_cast<CanvasRenderer*>(comboCanvases_->itemData(index).value<void*>());
    if (canvasRenderer_ != 0) {
        canvasSize_ = canvasRenderer_->getCanvasSize();

        int width = canvasSize_.x;
        if (width % 2 == 1)
        width -= 1;

        int height = canvasSize_.y;
        if (height % 2 == 1)
            height -= 1;

        spinWidth_->setValue(width);
        spinHeight_->setValue(height);
    }

}

void AnimationExportWidget::createWidgets(){

    mainLayout_ = new QVBoxLayout();
    renderBox_ = createAnimationRenderBox(this);
    mainLayout_->addWidget(renderBox_);

    setLayout(mainLayout_);
    setWidgetState();
}

void AnimationExportWidget::closeEvent(QCloseEvent* e) {
    currentFrame_ = static_cast<int>(duration_ / fpsFactor_ + 1);
    endRendering();
    QDialog::closeEvent(e);
}

void AnimationExportWidget::createConnections() {
    connect(comboCanvases_, SIGNAL(currentIndexChanged(int)), this, SLOT(controlledCanvasChanged(int)));
    connect(saveAsFrameSequenceButton_, SIGNAL(clicked()), this, SLOT(recordAnimationFrameSeq()));
    connect(saveAsVideoButton_, SIGNAL(clicked()), this, SLOT(recordAnimationVideo()));
    connect(videoSetupButton_, SIGNAL(clicked()), this, SLOT(videoSetup()));
}

void AnimationExportWidget::setWidgetState() {
    comboCanvases_->setEnabled(true);
}

void AnimationExportWidget::refreshComboBoxes() {
    if (comboCanvases_ != 0) {
        comboCanvases_->clear();
        for (std::vector<CanvasRenderer*>::const_iterator it = allCanvasRenderer_.begin(); it != allCanvasRenderer_.end(); ++it)
            comboCanvases_->addItem((*it)->getID().c_str(),qVariantFromValue(static_cast<void*>(*it)));
    }
}

QGroupBox* AnimationExportWidget::createAnimationRenderBox(QWidget* parent) {
    QGroupBox* recordBox = new QGroupBox(tr("Recording"), parent);
    QVBoxLayout* layout = new QVBoxLayout();
    QHBoxLayout* rowLayout = new QHBoxLayout();

    comboCanvases_ = new QComboBox(recordBox);
    rowLayout->addWidget(new QLabel(tr("Canvas:"), recordBox));
    rowLayout->addWidget(comboCanvases_);
    layout->addLayout(rowLayout);

    // render meta settings
    rowLayout = new QHBoxLayout();
    rowLayout->addWidget(new QLabel(tr("FPS:")));
    spinRecordingFPS_ = new QSpinBox(recordBox);
    spinRecordingFPS_->setRange(1, 120);
    spinRecordingFPS_->setValue(25);
    rowLayout->addWidget(spinRecordingFPS_);
    rowLayout->addStretch();

    rowLayout->addWidget(new QLabel(("     Video Dimensions:")), 1, 0);
    spinWidth_ = new QSpinBox(recordBox);
    spinHeight_ = new QSpinBox(recordBox);

    spinWidth_->setRange(64, 2048);
    spinWidth_->setSingleStep(4);
    spinWidth_->setValue(512);
    spinWidth_->setAccelerated(true);

    spinHeight_->setRange(64, 2048);
    spinHeight_->setSingleStep(4);
    spinHeight_->setValue(512);
    spinHeight_->setAccelerated(true);

    rowLayout->addWidget(spinWidth_);
    rowLayout->addWidget(new QLabel(" x ", 0, 0));
    rowLayout->addWidget(spinHeight_);
    layout->addLayout(rowLayout);

    QTime temp(0, 0, 0, 0);
    QTime duration = temp.addSecs((int)duration_/30);
    QHBoxLayout* row2Layout = new QHBoxLayout();
    spinStartTime_ = new QTimeEdit(QTime(0, 0 , 0, 0));
    spinEndTime_ = new QTimeEdit(duration);
    spinStartTime_->setDisplayFormat("mm:ss");
    spinEndTime_->setDisplayFormat("mm:ss");
    row2Layout->addWidget(new QLabel("Start Time: "));
    row2Layout->addWidget(spinStartTime_);

    spinStartTime_->setMaximumTime(duration);
    spinEndTime_->setMaximumTime(duration);

    row2Layout->addWidget(new QLabel("End Time: "));
    row2Layout->addWidget(spinEndTime_);
    layout->addLayout(row2Layout);

    // render 2 vid or frameSeq
    rowLayout = new QHBoxLayout();
    saveAsFrameSequenceButton_ = new QPushButton(tr("Save as frame sequence"));
    saveAsVideoButton_ = new QPushButton(tr("Save as video file"));
    videoSetupButton_ = new QPushButton(tr("Setup video..."));
#ifndef VRN_MODULE_FFMPEG
    saveAsVideoButton_->setVisible(false);
    videoSetupButton_->setVisible(false);
#endif

    rowLayout->addWidget(saveAsFrameSequenceButton_);
    rowLayout->addWidget(saveAsVideoButton_);
    rowLayout->addWidget(videoSetupButton_);
    layout->addLayout(rowLayout);

    recordBox->setLayout(layout);
    return recordBox;
}

void AnimationExportWidget::recordAnimationFrameSeq() {
    recordAnimation(false);
}

void AnimationExportWidget::recordAnimationVideo() {
    recordAnimation(true);
}

void AnimationExportWidget::recordAnimation(bool recordVideo) {
    renderingVideo_ = recordVideo;
    QStringList fileList;
    if (renderingVideo_ == false) {
        renderState_= Snapshot;
        QString s = QFileDialog::getExistingDirectory(this, tr("Select Output Directory"), VoreenApplication::app()->getUserDataPath().c_str());

        if (!s.isEmpty())
            fileList.push_back(s);
    }
#ifdef VRN_MODULE_FFMPEG
    else {
        renderState_= Recording;
        std::vector<std::string> formats = VideoEncoder::getSupportedFormatsByFileEnding();
        std::vector<std::string> descs = VideoEncoder::getSupportedFormatDescriptions();
        QStringList lstFormats;
        for (size_t i = 0; i < formats.size(); ++i) {
            formats[i] = (descs[i] + " (*." + formats[i] + ")");
            lstFormats.append(formats[i].c_str());
        }

        QFileDialog dialog(this);
        dialog.setDefaultSuffix("avi");
        dialog.setWindowTitle(tr("Export As Video File"));
        dialog.setDirectory(VoreenApplication::app()->getUserDataPath().c_str());
        #if QT_VERSION >= 0x040400
        dialog.setNameFilters(lstFormats);
        #endif
        dialog.setAcceptMode(QFileDialog::AcceptSave);

        if (dialog.exec())
            fileList = dialog.selectedFiles();
    }
#endif // VRN_MODULE_FFMPEG
    if (fileList.size() > 0) {
        recordPathName_ = fileList.first().toStdString();
        startRendering();
    }
}

void AnimationExportWidget::startRendering() {
    // gather some useful params
    fps_ = spinRecordingFPS_->value();
    startframe_= fps_ * (spinStartTime_->time().hour() * 3600 + spinStartTime_->time().minute()*60+ spinStartTime_->time().second());
    endframe_ = fps_ * (spinEndTime_->time().hour() * 3600 + spinEndTime_->time().minute()*60+ spinEndTime_->time().second());
    currentFrame_ = startframe_;
    duration_= endframe_;

    // set accurate recording dimension
    canvasSize_ = canvasRenderer_->getCanvasSize();
    dynamic_cast<tgt::QtCanvas*>(canvasRenderer_->getCanvas())->resize(spinWidth_->value(), spinHeight_->value());
    qApp->processEvents();

#ifdef VRN_MODULE_FFMPEG
    if (renderState_== Recording) {
        tgt::Texture* texture_ = canvasRenderer_->getImageColorTexture();
        try {
            ffmpegEncoder_.startVideoEncoding(recordPathName_.c_str(), fps_, spinWidth_->value(), spinHeight_->value(),
                                              texture_->getGLFormat(), texture_->getGLDataType());
        }
        catch (tgt::Exception& e) {
            QMessageBox::critical(this, tr("Video Export Failed"),
                                  tr("Failed to initialize video export:\n%1").arg(e.what()));
            return;
        }
    }
#endif // VRN_MODULE_FFMPEG

    VoreenQualityMode::getRef().requestQualityModeForced(VoreenQualityMode::RQ_HIGH, this);

    setWidgetState();
    dynamic_cast<tgt::QtCanvas*>(canvasRenderer_->getCanvas())->resize(spinWidth_->value(), spinHeight_->value());
    QProgressDialog progress(tr("Exporting Animation..."), tr("Abort"), 0, (int)(duration_ /fpsFactor_), this);
    progress.setWindowTitle(tr("Exporting"));
    if (fps_ > 0) {
        for(int i = 0; i < duration_ /fpsFactor_; ++i) {
            renderingStep();
            progress.setValue(i);
            qApp->processEvents();

            if (currentFrame_ > duration_ / fpsFactor_ || progress.wasCanceled()){
                break;
            }
        }
        endRendering();
    }
}

// rotate the camera about an arbitrary axis and angle
void AnimationExportWidget::rotateView(CameraProperty* camProp, float angle, const tgt::vec3& axis, const tgt::vec3& camPos, const tgt::vec3& camLook) {
    float cosA = cos(angle);
    float sinA = sin(angle);
    float ux2 = axis[0]*axis[0];
    float uy2 = axis[1]*axis[1];
    float uz2 = axis[2]*axis[2];
    float ux_uy = axis[0]*axis[1];
    float ux_uz = axis[0]*axis[2];
    float uy_uz = axis[1]*axis[2];
    tgt::mat3 mat(cosA + ux2*(1-cosA), ux_uy*(1-cosA)-axis[2]*sinA, ux_uz*(1-cosA)+axis[1]*sinA, ux_uy*(1-cosA)+axis[2]*sinA, cosA+uy2*(1-cosA), uy_uz*(1-cosA)-axis[0]*sinA, ux_uz*(1-cosA)-axis[1]*sinA, uy_uz*(1-cosA)+axis[0]*sinA, cosA+uz2*(1-cosA));
    tgt::vec3 look = mat*camLook;

    // set new focus-point
    camProp->setFocus(camPos + look);
}

void AnimationExportWidget::renderingStep(){
    if (animation_ == 0)
        return;
    tgtAssert(canvasRenderer_, "No canvas renderer");

    if (currentFrame_ > duration_ / fpsFactor_){
        endRendering();
    }
    else {
        animation_->renderAt((float)currentFrame_*fpsFactor_/ fps_);
        #ifdef VRN_MODULE_FFMPEG
        if (renderState_== Recording) {
            if (canvasRenderer_->getCanvasSize() != tgt::ivec2(spinWidth_->value(), spinHeight_->value())) {
                tgt::QtCanvas* canvas = dynamic_cast<tgt::QtCanvas*>(canvasRenderer_->getCanvas());
                canvas->resize(spinWidth_->value(), spinHeight_->value());
                canvas->repaint();
            }

            tgt::Texture* texture = canvasRenderer_->getImageColorTexture();
            if (texture && texture->getDimensions().xy() == tgt::ivec2(spinWidth_->value(), spinHeight_->value())) {
                texture->downloadTexture();
                ffmpegEncoder_.nextFrame(texture->getCpuTextureData());
            }
            else {
                LERRORC("voreenqt.AnimationExportWidget",
                        "Frame texture could not be downloaded or dimensions do not match");
            }
        } else {
        #endif
            // render frame to file
            char fn[1024];

            sprintf(fn, "%s%05d%s", std::string(recordPathName_ + "/frame").c_str(), currentFrame_, ".png");
            try {
                canvasRenderer_->renderToImage(fn, tgt::ivec2(spinWidth_->value(), spinHeight_->value()));
            }
            catch (const VoreenException& e) {
                LERRORC("voreen.AnimationExportWidget", "Failed to write image: " << std::string(fn) << " with error: " << e.what());
            }
        #ifdef VRN_MODULE_FFMPEG
        }
        #endif
        ++currentFrame_;
    }
}

void AnimationExportWidget::endRendering() {
    VoreenQualityMode::getRef().requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
    renderState_= Inactive;
#ifdef VRN_MODULE_FFMPEG
    if(renderingVideo_){
        ffmpegEncoder_.stopVideoEncoding();
        renderingVideo_ = false;
    }
#endif
    dynamic_cast<tgt::QtCanvas*>(canvasRenderer_->getCanvas())->resize(canvasSize_.x, canvasSize_.y);
}

void AnimationExportWidget::videoSetup() {
#ifdef VRN_MODULE_FFMPEG
    int curPreset = ffmpegEncoder_.getPreset();
    int curBitrate = ffmpegEncoder_.getBitrate();
    QDialog* dialog = createVideoSetupDialog(this, curPreset, curBitrate);
    if (dialog->exec() == QDialog::Accepted)
        ffmpegEncoder_.setup(preset_->currentIndex(), bitrate_->value() * 1024);

    delete dialog;
#endif
}

#ifdef VRN_MODULE_FFMPEG

QDialog* AnimationExportWidget::createVideoSetupDialog(QWidget* parent, int curPreset, int curBitrate) {
    QDialog* dialog = new QDialog(parent);
    dialog->setWindowTitle(tr("Setup Video"));
    QVBoxLayout* layout = new QVBoxLayout();
    dialog->setLayout(layout);

    // Preset
    QHBoxLayout* rowLayout = new QHBoxLayout();
    preset_ = new QComboBox(dialog);
    rowLayout->addWidget(new QLabel(tr("Preset:")));

    rowLayout->addWidget(preset_);
    const char** ccPairNames = VideoEncoder::getContainerCodecPairNames();
    for (int i = VideoEncoder::GUESS; i < VideoEncoder::LAST; ++i)
        preset_->addItem(ccPairNames[i]);
    preset_->setCurrentIndex(curPreset);
    layout->addLayout(rowLayout);

    // Bitrate
    rowLayout = new QHBoxLayout();
    bitrate_ = new QSpinBox(dialog);
    rowLayout->addWidget(new QLabel(tr("Bitrate (kbit/s):")));
    rowLayout->addWidget(bitrate_);
    bitrate_->setMinimum(400);
    bitrate_->setMaximum(4000);
    bitrate_->setSingleStep(1);
    bitrate_->setValue(curBitrate / 1024);
    layout->addLayout(rowLayout);

    rowLayout = new QHBoxLayout();
    QPushButton* ok = new QPushButton(QObject::tr("OK"));
    QPushButton* cancel = new QPushButton(QObject::tr("Cancel"));

    connect(ok, SIGNAL(clicked()), dialog, SLOT(accept()));
    connect(cancel, SIGNAL(clicked()), dialog, SLOT(reject()));

    rowLayout->addWidget(ok);
    rowLayout->addWidget(cancel);
    layout->addLayout(rowLayout);

    return dialog;
}

#endif // VRN_MODULE_FFMPEG

} // namespace voreen
