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

#ifndef VRN_CANVASRENDERERWIDGET_H
#define VRN_CANVASRENDERERWIDGET_H

#include "modules/core/processors/output/canvasrenderer.h"

#include "voreen/qt/widgets/processor/qprocessorwidget.h"

#include "tgt/qt/qtcanvas.h"
#include "voreen/qt/voreenqtapi.h"

namespace voreen {

class NetworkEvaluator;

/**
 * Widget of the CanvasRenderer.
 * It handles the transformation of touch events.
 */
class VRN_QT_API CanvasRendererWidget : public QProcessorWidget {
Q_OBJECT
public:
    /** Constructor */
    CanvasRendererWidget(QWidget* parent, CanvasRenderer* canvasRenderer);
    /** Destructor */
    virtual ~CanvasRendererWidget();
    /** Initializes the OpenGL */
    virtual void initialize();
    /** Updates all relevant properties and meta data */
    virtual void updateFromProcessor();
protected:
    /** Handling of touch events. */
    bool event(QEvent *event);
    /** Widget can be enlarged to fullscreen by pressing F11. */
    void keyPressEvent(QKeyEvent*);
    /** Emits a sizeChanged signal. */
    void resizeEvent(QResizeEvent*);

    tgt::QtCanvas* canvasWidget_;                   ///< canvas also used by stereocanvasrendererwidget

signals:
    void canvasRendererWidget_sizeChanged_Signal(); ///< emitted on resizeEvent
};

} // namespace voreen

#endif
