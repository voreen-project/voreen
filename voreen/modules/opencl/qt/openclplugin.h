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

#ifndef VRN_OPENCLPLUGIN_H
#define VRN_OPENCLPLUGIN_H

#include "voreen/core/voreenapplication.h"
#include "voreen/qt/widgets/codeedit.h"
#include "modules/opencl/qt/openclhighlighter.h"
#include "modules/opencl/properties/openclproperty.h"

#include <QToolButton>
#include <QUrl>
#include <QWidget>
#include <QSpinBox>

class QTextEdit;

namespace voreen {

class OpenCLPlugin : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor
     *
     * @param prop the shader property that belongs to this plugin
     * @param parent the parent widget
     */
    OpenCLPlugin(OpenCLProperty* prop, QWidget* parent = 0);

    /**
     * Destructor
     */
    ~OpenCLPlugin();

    /**
     * Creates all necessary widgets.
     */
    void createWidgets();

    /**
     * Creates all necessary connections.
     */
    void createConnections();
    void updateFromProperty();

signals:
    void modified();

public slots:
    void undoProgram();
    void fullUndoProgram();
    void openProgram();
    void saveProgram();
    void setProperty();
    void changeFontSize();

private:
    const QString getOpenFileName(QString filter);
    const QString getSaveFileName(QStringList filters);

    OpenCLProperty* property_;               ///< shader property that belongs to this plugin
    QToolButton* undoBt_;
    QToolButton* fullUndoBt_;
    QToolButton* openBt_;
    QToolButton* saveBt_;
    QToolButton* updateBt_;
    QSpinBox* fontSizeBox_;
    CodeEdit* codeEdit_;
    QTextEdit* compilerLogWidget_;
    OpenCLHighlighter* highlighter_;
};

} // namespace voreen

#endif // VRN_OPENCLPLUGIN_H
