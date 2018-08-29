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

#ifndef VRN_TRANSFUNCPLUGINXXX_H
#define VRN_TRANSFUNCPLUGINXXX_H

#include <QWidget>

class QTabWidget;
class QShowEvent;

namespace voreen {

class TransFuncEditorBase;
class TransFuncProperty;

/**
 * Container class for transfer function editors. Actually 4 editors for intensity
 * transfer functions are possible. They are arranged in a QTabWidget.
 */
class TransFuncPluginXXX : public QWidget {
    Q_OBJECT
public:

    /**
     * Constructor
     *
     * @param prop the transfer function property that belongs to this plugin
     * @param parent the parent widget
     * @param orientation should the editors layouted vertically or horizontally?
     */
    TransFuncPluginXXX(TransFuncProperty* prop, QWidget* parent = 0,
                    Qt::Orientation orientation = Qt::Horizontal);

    /**
     * Destructor
     */
    ~TransFuncPluginXXX();

    /**
     * Creates the editors
     *
     * @param parent the parent widget
     * @param orientation should the editors layouted vertically or horizontally?
     */
    void createEditors(QWidget* parent, Qt::Orientation orientation);

    /**
     * Creates the widgets for all editors and puts them into a QTabWidget.
     */
    void createWidgets();

    /**
     * Creates the connections for all editors and the QTabWidget.
     * Additionally the first editor is set active by calling editorChanged()
     */
    void createConnections();

    /**
     * This method is called by the transfer function property when the rendered volume changes.
     * It calls the volumeChanged() method of the current active editor.
     */
    void updateFromProperty();

    void disconnect();

public slots:
    /**
     * Slot that is called when the user selects another editor in the tab widget.
     * It will switch to the new editor when the user agrees with the reset of the
     * transfer function. Otherwise the old editor remains.
     *
     * @param index index of the selected tab
     */
    void editorChanged(int index);

signals:
    /**
     * This signal is emitted whenever the transfer function has changed. It is used to force a
     * repaint of the volume rendering
     */
    void transferFunctionChanged();

protected:
    /**
     * Called by Qt when the widget is shown. Calls update() on the current active TransFuncEditor.
     *
     * @param event the showevent
     */
    void showEvent(QShowEvent* event);

private:
    std::vector<TransFuncEditorBase*> editors_; ///< vector with all editors. The user can choose between all.
    TransFuncProperty* property_;               ///< transfer function property that belongs to this plugin
    QTabWidget* tabWidget_;                 ///< QTabWidget with all created editors
    int oldIndex_;                          ///< index of the last active tab
};

} // namespace voreen

#endif // VRN_TRANSFUNCPLUGIN_H
