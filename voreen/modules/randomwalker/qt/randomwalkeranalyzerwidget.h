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

#ifndef VRN_RANDOMWALKERANALYZERWIDGET_H
#define VRN_RANDOMWALKERANALYZERWIDGET_H

#include "voreen/qt/widgets/processor/qprocessorwidget.h"
#include "../processors/randomwalkeranalyzer.h"

#include <QWidget>
#include <QPushButton>
#include <QToolButton>
#include <QSpinBox>
#include <QLabel>
#include <QFileDialog>
#include <QLineEdit>
#include <QSlider>

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsEllipseItem>
#include <QTableWidget>

#include "tgt/event/eventlistener.h"

namespace voreen {

class RandomWalkerAnalyzerWidget : public QProcessorWidget {
    Q_OBJECT
public:
    RandomWalkerAnalyzerWidget(RandomWalkerAnalyzer* processor, QWidget* parent);
    ~RandomWalkerAnalyzerWidget();

    virtual void updateFromProcessor();

    virtual void processorNameChanged();

protected slots:
    void computeRandomWalker();
    void zoomIn();
    void resetZoom();
    void tableSelectionChanged();

private:

    void updateTableWidget();

    QTableWidget* tableUncertaintyRegions_;
    QPushButton* computeButton_;
    QPushButton* zoomInButton_;
    QPushButton* resetZoomButton_;

    bool blockProcessorUpdates_;

    int curFocusRegion_;

    static const std::string loggerCat_;

};

} // namespace voreen

#endif // VRN_RANDOMWALKERANALYZERWIDGET_H
