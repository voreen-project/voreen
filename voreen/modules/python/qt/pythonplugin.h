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

#ifndef VRN_PYTHONPLUGIN_H
#define VRN_PYTHONPLUGIN_H

#include <QWidget>
#include "../core/pythonoutputlistener.h"

class QTextEdit;
class QToolButton;
class CodeEdit;

namespace voreen {

class PythonHighlighter;
class PythonScript;
class PythonProperty;

class PythonPlugin : public QWidget, public PythonOutputListener {
    Q_OBJECT
public:

    /**
     * Constructor
     *
     * @param script the python property that belongs to this plugin.
     * @param parent the parent widget
     */
    PythonPlugin(PythonProperty* property = nullptr, QWidget* parent = nullptr);
    ~PythonPlugin();

    void clearScript();
    void updateFromProperty();

    /// Output listener functions.
    virtual void pyStdout(const std::string& out, const std::string& id);
    virtual void pyStderr(const std::string& err, const std::string& id);

public slots:
    void runScript();
    void newScript();

signals:
    void modified();

private:
    const QString selectOpenFileName(QString filter);
    const QString selectSaveFileName(QStringList filters);
    bool saveScriptInternal(QString filename, QString source);
    void tryLoadScript(QString fileName);

private slots:
    void openScript();
    void reloadScript();
    void saveScript();
    void saveScriptAs();

    void increaseFontSize();
    void decreaseFontSize();
    void updateFont();

    void updateGuiState();

private:
    PythonScript* script_;
    PythonProperty* property_;

    QToolButton* runBt_;
    QToolButton* newBt_;
    QToolButton* openBt_;
    QToolButton* reloadBt_;
    QToolButton* saveBt_;
    QToolButton* saveAsBt_;

    CodeEdit* codeEdit_;
    QTextEdit* compilerLogWidget_;
    PythonHighlighter* highlighter_;

    int fontSize_;
    bool scriptOwner_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PYTHONEDITOR_H
