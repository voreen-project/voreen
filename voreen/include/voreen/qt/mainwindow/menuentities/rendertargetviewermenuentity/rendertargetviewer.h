/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_RENDERTARGETVIEWER_H
#define VRN_RENDERTARGETVIEWER_H

#include <vector>

#include <QAction>
//#include <QtOpenGL/QGLWidget>

#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/workspace.h"

#include "voreen/qt/voreenqtapi.h"

namespace tgt {
    class QtCanvas;
}

namespace voreen {

class RenderTarget;

/**
 * Show currently used RenderTargets for debugging purposes.
 */
class VRN_QT_API RenderTargetViewer : public QWidget, public NetworkEvaluatorObserver, public UsesWorkspace {
Q_OBJECT
    friend class RenderTargetViewerPainter;
    friend class RenderTargetViewerMenuEntity;
public:
    enum ShowType {
        R                   = 0x0001,
        G                   = 0x0002,
        B                   = 0x0004,
        A                   = 0x0008,
        Color               = 0x0010,
        Alpha               = 0x0020,
        Depth               = 0x0040,
        CheckerboardPattern = 0x0080,
        BackgroundWhite     = 0x0100,
        BackgroundBlack     = 0x0200,
        H                   = 0x0400,
        S                   = 0x0800,
        V                   = 0x1000
    };

    RenderTargetViewer();

    ~RenderTargetViewer();

    void setEvaluator(NetworkEvaluator* evaluator);

    /** @override NetworkEvaluatorObserver */
    void afterNetworkProcess();
    /** @override UsesWorkspace */
    void setWorkspace(Workspace* workspace);

    void initialize();
    void deinitialize();

    //--------------------------
    //  helper functions
    //--------------------------
protected:
    /**
     * Collects all render ports from the current network.
     */
    std::vector<RenderPort*> collectRenderPorts();

    void takeScreenshot(bool withOverlay);

public slots:
    /**
     * Connected to NetworkEditor.selectedProcessor
     * Connection done in RenderTargetViewerMenuEntity.
     */
    void processorsSelected(const QList<Processor*>& processors);

    //--------------------------
    //  event handling
    //--------------------------
protected:
    void mousePressEvent(QMouseEvent* e);
    void mouseReleaseEvent(QMouseEvent* e);
    void mouseMoveEvent(QMouseEvent* e);
    void wheelEvent(QWheelEvent* e);
    void keyPressEvent(QKeyEvent* e);
    void closeEvent(QCloseEvent* e);

    //--------------------------
    //  event helper
    //--------------------------
    void updateMousePosition(QMouseEvent* e);
    void updateSelected();

    //--------------------------
    //  member
    //--------------------------
protected:
    NetworkEvaluator* evaluator_; ///< NetworkEvaluator to access render ports
    tgt::QtCanvas* canvas_;       ///< Canvas used to paint the render ports

    QList<Processor*> selectedProcessors_; ///< list of selected processors of the network
    int selectedRenderPortIndex_;          ///< -1 if no port is selected
    bool maximizeOnePort_;                 ///< true, if only one port is visualized
    std::map<int, unsigned int> showType_; ///< map storing all selected options associated with a render port

    //helper used in the painter
    int mouseX_;            ///< current mouse x position
    int mouseY_;            ///< current mouse y position
    bool mouseIsInside_;

    //zoom helper used in the painter
    float zoomScale_;
    float zoomTranslateX_;
    float zoomTranslateY_;
    float zoomMouseX_;
    float zoomMouseY_;
    float zoomOffsetX_;
    float zoomOffsetY_;

    static const std::string loggerCat_;
protected:
    //Qt actions. TODO: needed?
    QMenu* contextMenuMEN_;
    QActionGroup* typeToShowACG_;
    QActionGroup* backgroundToShowACG_;
    QAction* backgroundBlackACT_;
    QAction* backgroundWhiteACT_;
    QAction* backgroundCheckerboardPatternACT_;
    QAction* colorRACT_;
    QAction* colorGACT_;
    QAction* colorBACT_;
    QAction* colorAACT_;
    QAction* colorRGBAACT_;
    QAction* alphaOnlyACT_;
    QAction* depthOnlyACT_;
    QAction* hOnlyACT_;
    QAction* sOnlyACT_;
    QAction* vOnlyACT_;
    QAction* keepAspectRatioACT_;
    QAction* showInfosACT_;
    QAction* showInfosDetailsACT_;
    QAction* saveScreenshotACT_;
    QAction* saveScreenshotWithOverlayACT_;
    QAction* filterPortyBySelectedProcessorsACT_;
};

} // namespace voreen

#endif //VRN_RENDERTARGETVIEWER_H
