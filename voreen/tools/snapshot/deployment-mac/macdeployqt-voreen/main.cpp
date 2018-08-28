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

#include "shared/shared.h"

int main(int argc, char **argv)
{
    QString appBundlePath;
    if (argc > 1)
        appBundlePath = QString::fromLocal8Bit(argv[1]);

    if (argc < 2 || appBundlePath.startsWith("-")) {
        qDebug() << "Usage: macdeployqt-voreen app-bundle [-no-frameworks] [-no-plugins] [-no-dylibs] [-dmg]";
        qDebug() << "";
        qDebug() << "macdeployqt-voreen creates a self-contained application bundle that";
        qDebug() << "contains the Qt frameworks, plugins and dylibs used by the application.";
        qDebug() << "";
        qDebug() << "macdeployqt-voreen is based on the Qt macdeployqt tool,";
        qDebug() << "but additionally deploys dylibs used by the application.";
        qDebug() << "";
        qDebug() << "Qt-Frameworks in use are copied into the bundle, unless \"-no-frameworks\" is specified.";
        qDebug() << "Plugins related to a framework are copied in with the";
        qDebug() << "framework. The accessibilty, image formats, and text codec";
        qDebug() << "plugins are always copied, unless \"-no-plugins\" is specified. ";
        qDebug() << "Additional dylibs are copied, unless \"-no-dylibs\" is specified.";
        qDebug() << "";
        qDebug() << "See the \"Deploying an Application on Qt/Mac\" typic in the";
        qDebug() << "documentation for more information about deployment on Mac OS X.";

        return 0;
    }
    
    if (appBundlePath.endsWith("/"))
        appBundlePath.chop(1);
    
    bool frameworks = true;
    bool plugins = true;
    bool dylibs = true;
    bool dmg = false;

    for (int i = 1; i < argc; ++i) {
        QByteArray argument = QByteArray(argv[i]);
        if (argument == QByteArray("-no-frameworks"))
            frameworks = false;
        if (argument == QByteArray("-no-plugins"))
            plugins = false;
        if (argument == QByteArray("-no-dylibs"))
            dylibs = false;
        if (argument == QByteArray("-dmg"))
            dmg = true;
    }

    if (frameworks) {
        qDebug() << "";
        qDebug() << "Deploying frameworks to" << appBundlePath;
        DeploymentInfo deploymentInfo  = deployQtFrameworks(appBundlePath);

        // cannot deploy plugins without deploying frameworks
        if (plugins) {
            if (deploymentInfo.qtPath.isEmpty())
                deploymentInfo.pluginPath = "/Developer/Applications/Qt/plugins"; // Assume binary package.
            else
                deploymentInfo.pluginPath = deploymentInfo.qtPath + "/plugins";

            qDebug() << "";
            qDebug() << "Deploying plugins from" << deploymentInfo.pluginPath;
            deployPlugins(appBundlePath, deploymentInfo);
            createQtConf(appBundlePath);
        }
    }

    if (dylibs) {
        qDebug() << "Deploying dylibs to" << appBundlePath;
        deployDylibs(appBundlePath);
    }

    if (dmg) {
        qDebug() << "";
        qDebug() << "Creating disk image (.dmg) for" << appBundlePath;
        createDiskImage(appBundlePath);
    }
}

