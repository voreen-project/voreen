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

#ifndef VRN_BODYPARTS3DWIDGET_H
#define VRN_BODYPARTS3DWIDGET_H

#include "voreen/qt/widgets/processor/qprocessorwidget.h"
#include "voreen/qt/voreenqtapi.h"

#include "modules/touchtable/processors/bodyparts3d/bodyparts3dsource.h"
#include "../datastructures/bodyparts3dmodel.h"

#include <QTreeView>
#include <QToolBar>
#include <QAction>
#include <QComboBox>
#include <QVBoxLayout>

namespace voreen{

class VRN_QT_API BodyParts3DWidget : public QProcessorWidget {
Q_OBJECT
public:
    BodyParts3DWidget(QWidget* parent, BodyParts3DSource* bodyParts3DSource);
    virtual ~BodyParts3DWidget();

    /**
     * This method updates the content of combo_, version_, path_, and creates a new model
     */
    virtual void updateFromProcessor();
protected slots:
    /**
     * sets the version according to the selection of combo_
     */
    void switchversion();

    /**
     * opens a filedialog to get a path and filename and passes it to the model to save the current Colorscheme
     */
    void save();

    /**
     * loads a custom or standard Color scheme for BodyParts3D.
     *
     * note that the path containint that file must be the same path as the folder of the geometry files.
     */
    void load();

    /**
     * loads the filenames, colors and descriptions for the checked BodyParts and passes them to the processor.
     */
    void loadBodyPartsData();

    /**
     * Opens a colorpicker after a click in the color column
     */
    void clicked(const QModelIndex& index);
private:
    QVBoxLayout* mainLayout_; ///< windowlayout
    QTreeView* view_; ///< view for the Atlas
    BodyParts3DModel* model_; ///< Atlas datastructure

    QComboBox* combo_; ///< dropdown menu for version selection

    QAction* saveAction_; ///< save Button
    QAction* loadAction_; ///< load Button
    QAction* loadBodyPartsDataAction_;
    QToolBar* toolBar_; ///< Toolbar for the Buttons

    BodyPartsVersion version_; ///< current BodyParts Version
    std::string path_; ///< current Path

    static const std::string loggerCat_;
};

}//namespace voreen

#endif // VRN_BODYPARTS3DWIDGET_H
