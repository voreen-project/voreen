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

#ifndef VRN_PYTHONEDITOR_H
#define VRN_PYTHONEDITOR_H

// include Python at very first
#include "modules/python/pythonmodule.h"

#include "modules/python/qt/pythonhighlighter.h"

#include "voreen/qt/widgets/codeedit.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h"

#include <QToolButton>
#include <QUrl>
#include <QWidget>
#include <QObject>

class QTextEdit;

namespace voreen {

class PythonScript;
class NetworkEvaluator;

class PythonEditor : public QObject, public VoreenQtMenuEntity, public PythonOutputListener {
    Q_OBJECT
public:
    PythonEditor();
    ~PythonEditor();

    virtual std::string getName() const { return "Python Scripting"; }
    virtual QIcon getIcon() const       { return QIcon(":/modules/python/python.png"); }

    void clearScript();

    /// Output listener functions.
    virtual void pyStdout(const std::string& out);
    virtual void pyStderr(const std::string& err);

signals:
    void modified();

protected:
    virtual QWidget* createWidget() const;

    virtual void initialize();
    virtual void deinitialize();

private:
    const QString selectOpenFileName(QString filter);
    const QString selectSaveFileName(QStringList filters);
    bool saveScriptInternal(QString filename, QString source);
    void tryLoadScript(QString fileName);

private slots:
    void runScript();
    void newScript();
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

    mutable QToolButton* runBt_;
    mutable QToolButton* newBt_;
    mutable QToolButton* openBt_;
    mutable QToolButton* reloadBt_;
    mutable QToolButton* saveBt_;
    mutable QToolButton* saveAsBt_;

    mutable QToolButton* increaseFontSizeBt_;
    mutable QToolButton* decreaseFontSizeBt_;

    mutable CodeEdit* codeEdit_;
    mutable QTextEdit* compilerLogWidget_;
    mutable PythonHighlighter* highlighter_;

    mutable QWidget* pythonWidget_;

    int fontSize_;
    bool scriptOwner_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PYTHONEDITOR_H
