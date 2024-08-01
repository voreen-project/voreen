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

#ifndef VRN_SHADERPLUGIN_H
#define VRN_SHADERPLUGIN_H

#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/qt/widgets/glslhighlighter.h"
#include "voreen/qt/widgets/codeedit.h"

#include <QToolButton>
#include <QUrl>
#include <QWidget>
#include <QSpinBox>
#include <QTabWidget>

class QTextEdit;

namespace voreen {

struct ShaderTab {
public:
    ShaderTab(CodeEdit* codeEdit, int sourceComponentIndex = NO_SOURCE_COMPONENT_INDEX);
    ShaderTab(ShaderTab&&);

    CodeEdit* codeEdit_;
    std::unique_ptr<GLSLHighlighter> highlighter_;
    int sourceComponentIndex_; /// Index of the source of this tab in the ShaderSource object of the property

    bool showsImported() const;

    static const int NO_SOURCE_COMPONENT_INDEX;
};

class ShaderPlugin : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor
     *
     * @param prop the shader property that belongs to this plugin
     * @param parent the parent widget
     */
    ShaderPlugin(ShaderProperty* prop, QWidget* parent = 0);

    /**
     * Destructor
     */
    ~ShaderPlugin();

    /**
     * Creates all necessary widgets.
     */
    void createWidgets();

    /**
     * Creates all necessary connections.
     */
    void createConnections();
    void updateFromProperty();

    void openImportTab(const std::string& filename);
    void openSourceComponentTab(int sourceComponentIndex);

signals:
    void modified();

public slots:
    void undoShader();
    void fullUndoShader();
    void openShader();
    void exportShader();
    void saveShader();
    void saveAllShaders();
    void changeFontSize();
    void urlClicked(const QUrl& link);
    void closeTab(int index);
    void adaptButtonsToTab(int index);

private:
    const QString getOpenFileName(QString filter);
    const QString getSaveFileName(QStringList filters);
    QString getFilter(tgt::ShaderObject::ShaderType type) const;
    ShaderTab& getCurrentShaderTab();
    void openTab(const std::string& filename, ShaderTab&& tab);
    void updateTabTexts();


    ShaderProperty* property_;               ///< shader property that belongs to this plugin
    QToolButton* undoBt_;
    QToolButton* fullUndoBt_;
    QToolButton* openBt_;
    QToolButton* exportBt_;
    QToolButton* updateBt_;
    QToolButton* updateAllBt_;
    QSpinBox* fontSizeBox_;

    QTabWidget* tabWidget_;
    QTextEdit* compilerLogWidget_;

    std::list<ShaderTab> shaderTabs_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_SHADERPLUGIN_H
