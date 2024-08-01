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

#ifdef WIN32
#include <windows.h>
#endif

#include <cstdlib>

#include "tgt/tgt_gl.h"

#include <GL/freeglut.h>

#include "tgt/shadermanager.h"
#include "tgt/logmanager.h"
#include "tgt/glcontextmanager.h"
#include "tgt/glut/glutcanvas.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/workspace.h"
#include "voreen/core/network/processornetwork.h"

#include "modules/core/processors/output/canvasrenderer.h"

const std::string description =
    "This is a small program that demonstrates a little of what the Voreen core library can do: \n"
    "We load the standard.vws workspace, which shows a volume rendering of the nucleon dataset \n"
    "that can be rotated and zoomed using the mouse.\n\n"
    "This is the GLUT-version of this sample, there are others in the simple/-folder, like a Qt-version.";

const std::string defaultWorkspace = "../voreenve/workspaces/rendering-raycasting.vws";

using namespace voreen;

//---------------------------------------------------------------------------

tgt::GLUTCanvas* canvas = 0;
VoreenApplication* app = 0;

NetworkEvaluator* networkEvaluator = 0;
ProcessorNetwork* network = 0;
VoreenPainter* painter = 0;

//---------------------------------------------------------------------------

void initialize() {
    Workspace* workspace = new Workspace();
    try {
        workspace->load(VoreenApplication::app()->getCoreResourcePath(defaultWorkspace));
    }
    catch (SerializationException& e) {
        LERRORC("simple-glut.initialize", "Failed to load standard workspace: " << e.what());
        exit(EXIT_FAILURE);
    }
    // initialize the network evaluator and retrieve CanvasRenderer processors from the loaded network
    networkEvaluator = new NetworkEvaluator();
    network = workspace->getProcessorNetwork();
    std::vector<CanvasRenderer*> canvasRenderer = network->getProcessorsByType<CanvasRenderer>();
    if (canvasRenderer.empty()) {
        LERRORC("simple-glut.initialize", "Loaded standard workspace does not contain a CanvasRenderer");
        exit(EXIT_FAILURE);
    }

    // init painter and connect it to the canvas
    painter = new VoreenPainter(canvas, networkEvaluator, canvasRenderer[0]);
    canvasRenderer[0]->setCanvas(canvas);
    // set the network to the network evaluator
    networkEvaluator->setProcessorNetwork(network);
}

void finalize() {
    delete painter;
    painter = 0;
    delete network;
    network = 0;
    delete networkEvaluator;
    networkEvaluator = 0;

    if (app) {
        app->deinitializeGL();
        app->deinitialize();
    }
    delete app;
    app = 0;
}

void keyPressed(unsigned char key, int /*x*/, int /*y*/) {
    switch (key) {
        case '\033': // = ESC
            finalize();
            exit(0);
            break;
        case '1':
            glutReshapeWindow(128,128);
            break;
        case '2':
            glutReshapeWindow(256,256);
            break;
        case '3':
            glutReshapeWindow(512,512);
            break;
        case '4':
            glutReshapeWindow(1024, 1024);
            break;
    }
}

int main(int argc, char** argv) {
    VoreenApplication app("simple-GLUT", "Simple-GLUT", description, argc, argv, VoreenApplication::APP_ALL);
    app.initialize();

    glutInit(&argc, argv);

#ifndef VRN_OPENGL_COMPATIBILITY_PROFILE
    glutInitContextVersion(VRN_OPENGL_CONTEXT_MAJOR_VERSION, VRN_OPENGL_CONTEXT_MINOR_VERSION);
    //glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
    glutInitContextProfile(GLUT_CORE_PROFILE);
#endif

    // initialize canvas
    canvas = new tgt::GLUTCanvas("Voreen - The Volume Rendering Engine (Simple-GLUT)",
                                  tgt::ivec2(512, 512), tgt::GLCanvas::RGBADD);
    canvas->init();

    glutKeyboardFunc(keyPressed);

    app.initializeGL();
    initialize();

    glutMainLoop();

    delete canvas;
    return 0; // will never be reached for standard GLUT
}
