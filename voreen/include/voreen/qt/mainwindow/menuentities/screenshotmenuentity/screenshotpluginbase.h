/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_SCREENSHOTPLUGINBASE_H
#define VRN_SCREENSHOTPLUGINBASE_H

#include <QWidget>

#include "voreen/qt/voreenqtapi.h"
#include "voreen/core/network/workspace.h"

#include <vector>

#include <QStringList>

class QToolButton;
class QSpinBox;
class QComboBox;
class QLineEdit;

namespace voreen {

/**
 * This is the base class of screenshot plugins used in the ScreenshotMenuEntity.
 * @see CanvasRendererScreenshotPlugin
 * @see NetworkScreenshotPlugin
 *
 * !!! loadMetaDAta() must be called in subclass constructor !!!
 * !!! nativeSizeHasChangedSlot() must be connected in subclass constructor !!!
 */
class VRN_QT_API ScreenshotPluginBase : public QWidget {
Q_OBJECT
    friend class ScreenshotMenuEntity;
public:
    /** Constructor */
    ScreenshotPluginBase(QWidget* parent, Workspace* workspace);
    /** Destructor */
    ~ScreenshotPluginBase();

    //-----------------
    //  To Override
    //-----------------
protected:
    /** Function triggered by the screenshot button. */
    virtual void saveScreenshot(const QString& filename, int width, int height) = 0;

    //-----------------
    //  Layout
    //-----------------
protected:
    /** Creates the entire base layout */
    void createLayout();

    QComboBox* resolutionComboBox_;         ///< combo box containing all supported resolutions
        QStringList resolutions_;           ///< list containing all supported strings
    QSpinBox*  userDefinedWidthSpinBox_;    ///< width, if user-defined is selected
    QSpinBox*  userDefinedHeightSpinBox_;   ///< height, if user-defined is selected
    QLineEdit* prefixLineEdit_;             ///< Used to set the save prefix
    QToolButton* screenshotButton_;         ///< button triggering the screenshot
        QString path_;                      ///< last used file

    //-----------------
    //  Callbacks
    //-----------------
public slots:
    virtual void nativeSizeHasChangedSlot() = 0;
    void resolutionComboBoxChangedSlot(int index);
    void userDefinedWidthSpinBoxChangedSlot(int value);
    void userDefinedHeightSpinBoxChangedSlot(int value);
    void prefixLineEditChangedSlot(QString prefix);
    void screenshotButtonPressedSlot();
    //----------------
    //  Meta Handling
    //----------------
protected:
    /** Updates the widget according to the workspace. */
    void updateWorkspace(Workspace* workspace);
    /** Loads the current meta information from the current workspace. */
    void loadMetaData();
    /** Stores the current meta information in the current workspace. */
    void saveMetaData();
    /** Helper used to remove old meta data. @see ScreenshotMenuEntity */
    void removeMetaData();

    Workspace* currentWorkspace_;       ///< currently used workspace
    std::string metaDataPrefix_;        ///< meta data prefix of the class
};

} // namespace voreen

#endif // VRN_SCREENSHOTPLUGINBASE_H
