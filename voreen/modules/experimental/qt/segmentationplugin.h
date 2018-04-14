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

#ifndef SEGMENTATIONPLUGIN_H
#define SEGMENTATIONPLUGIN_H

#include "../processors/regiongrowing.h"
#include "modules/base/processors/render/segmentationraycaster.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditor.h"
#include "tgt/event/eventlistener.h"

#include <QWidget>
#include <QCheckBox>

class QComboBox;
class QPushButton;
class QSpinBox;
class QDoubleSpinBox;
class QGroupBox;
class QComboBox;
class QLabel;

namespace voreen {

class NetworkEvaluator;
class ThresholdWidget;

class SegmentationPlugin : public QWidget, tgt::EventListener {
    Q_OBJECT
public:
    SegmentationPlugin(QWidget* parent, NetworkEvaluator* evaluator);
    ~SegmentationPlugin();

    void createWidgets();
    void createConnections();

    virtual void mousePressEvent(tgt::MouseEvent* e);

    bool usable(const std::vector<Processor*>& processors);

private slots:

    // update the intensity editor to the segmentation raycaster's state
    void updateIntensityEditor();

    // registers this plugin as listener at the segmentation raycaster's properties
    void registerAsListener();

    // is called when a relevant transfunc property was changed outside this plugin
    void transFuncChangedExternally();

    // callback for registered at segmentation raycaster's property
    void applySegmentationToggled();

    // callback gui checkbox
    void toggleApplySegmentation(bool);

    // gui callbacks
    void setCurrentSegment(int);
    void undoSegment();
    void clearSegment();
    void setSeed(bool checked);
    void saveSegmentation();
    void loadSegmentation();
    void clearSegmentation();

    void repaintCanvas();

private:

    // Propagates the locally saved region growing parameters to the RegionGrowing processor
    void propagateRegionGrowingParams();

    // Retrieve processors from evaluator
    SegmentationRaycaster* getSegmentationRaycaster(bool suppressWarning = false);
    RegionGrowingProcessor* getRegionGrowingProcessor(bool suppressWarning = false);

    NetworkEvaluator* evaluator_;
    TransFunc1DKeysPropertyEditor* intensityEditor_;

    QCheckBox* checkApplySegmentation_;
    QLabel* labelCurrentSegment_;
    QSpinBox* spinCurrentSegment_;

    QPushButton* markSeedButton_;
    QPushButton* undoButton_;
    QSpinBox* outputSegment_;
    QCheckBox* checkApplyThresholds_;
    QDoubleSpinBox* spinStrictness_;
    QComboBox* comboCostFunction_;
    QCheckBox* checkAdaptive_;
    QSpinBox* spinMaxSeedDist_;

    QPushButton* clearSegmentButton_;
    QPushButton* clearSegmentationButton_;
    QPushButton* saveSegmentationButton_;
    QPushButton* loadSegmentationButton_;

    QGroupBox* renderingBox_;
    QGroupBox* regionGrowingBox_;

};

} // namespace voreen

#endif
