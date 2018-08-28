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

#include <QApplication>
#include <QMainWindow>
#include <QMessageBox>

#include "tgt/init.h"
#include "tgt/logmanager.h"
#include "tgt/camera.h"
#include "tgt/shadermanager.h"
#include "tgt/glcontextmanager.h"
#include "tgt/qt/qtcanvas.h"

#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/qt/voreenapplicationqt.h"

#include "modules/core/processors/output/canvasrenderer.h"

#ifdef UNIX
#include <clocale>
#endif

using namespace voreen;

const std::string defaultWorkspace = "../voreenve/workspaces/raycasting-standard.vws";

const std::string description =
    "This is a small program that demonstrates a little of what the Voreen core library can do: \n"
    "We load the standard.vws workspace, which shows a volume rendering of the nucleon dataset \n"
    "that can be rotated and zoomed using the mouse.\n\n"
    "This is the Qt-version of this sample, there are others in the simple/-folder, like a GLUT-version.";

int main(int argc, char* argv[]) {

    // init both Qt and Voreen application
    VoreenApplicationQt::setupPrerequisites();
    QApplication myapp(argc, argv);
    VoreenApplicationQt vapp("simple-qt", "Simple-Qt", description, argc, argv,
        VoreenApplication::ApplicationFeatures(VoreenApplication::APP_ALL &~ VoreenApplication::APP_PROCESSOR_WIDGETS));

    // fixes problems with locale settings due to Qt (see http://doc.qt.io/qt-5/qcoreapplication.html#locale-settings)
#ifdef UNIX
    std::setlocale(LC_NUMERIC, "C");
#endif

    vapp.initialize();
    vapp.initializeGL();

    // create the mainwindow and assign a canvas to it as central widget
    QMainWindow* mainwindow = new QMainWindow();
    VoreenApplicationQt::qtApp()->setMainWindow(mainwindow);
    mainwindow->setWindowTitle("Voreen - The Volume Rendering Engine (Simple-Qt)");
    mainwindow->resize(512, 512);
    mainwindow->show();

    // load workspace from disc
    Workspace* workspace = new Workspace();
    try {
        workspace->load(VoreenApplication::app()->getCoreResourcePath(defaultWorkspace));
    }
    catch (SerializationException& e) {
        QMessageBox::critical(mainwindow, "Loading Workspace Failed", QString::fromStdString(e.what()));
        return EXIT_FAILURE;
    }

    // initialize the network evaluator and retrieve the CanvasRenderer processor from the loaded network
    NetworkEvaluator* networkEvaluator = new NetworkEvaluator();
    vapp.registerNetworkEvaluator(networkEvaluator);
    ProcessorNetwork* network = workspace->getProcessorNetwork();
    std::vector<CanvasRenderer*> canvasRenderer = network->getProcessorsByType<CanvasRenderer>();
    if (canvasRenderer.empty()) {
        QMessageBox::critical(mainwindow, "Invalid Workspace", "Loaded workspace does not contain a CanvasRenderer.");
        return EXIT_FAILURE;
    }

    // create the output canvas
    tgt::QtCanvas* canvas = new tgt::QtCanvas("Canvas", tgt::ivec2(256, 256), tgt::GLCanvas::RGBADD, 0);

    // init painter and connect it to canvas, evaluator and canvas renderer
    VoreenPainter* painter = new VoreenPainter(canvas, networkEvaluator, canvasRenderer[0]); //canvas takes ownership
    canvasRenderer[0]->setCanvas(canvas);

    // Add it to mainwindow.
    mainwindow->setCentralWidget(canvas);

    // Init has to be called after adding to layout!
    canvas->init();

    // pass the network to the network evaluator, which also initializes the processors
    networkEvaluator->setProcessorNetwork(network);

    // start the event process; the program runs as long as theres no exit-event
    myapp.exec();

    // we're done as soon as myapp.exec() returns, so we can delete everything
    delete workspace;
    delete networkEvaluator;
    delete mainwindow;

    vapp.deinitializeGL();
    vapp.deinitialize();

    return 0;
}
