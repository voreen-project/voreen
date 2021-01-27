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

#ifndef VRN_ABOUTBOXBASE_H
#define VRN_ABOUTBOXBASE_H

#include <QDialog>
#include <QStringList>

#include "voreen/qt/voreenqtapi.h"

namespace voreen {

    class VoreenQtMainWindow;

    /**
     * Base class of the about box.
     */
class VRN_QT_API AboutBoxBase : public QDialog {
    Q_OBJECT
public:
    /** Constructor
     *
     * @param mainWindow Pointer to the main window. Used to get application name and icon.
     * @param imagePath Filepath to the image shown in the about box
     */
    AboutBoxBase(VoreenQtMainWindow* mainWindow, QString imagePath = QString(""));

    /**
     * Calls "initialize" before exec
     * @override QDialog::exec()
     */
    virtual int exec();
private:
    //---------------------
    // basic layout
    //---------------------
    /**
     * Calls all init functions.
     * Is been called in exec() if not initialized.
     */
    void initialize();
    virtual void initAndLayoutItems();
    virtual void initSoftwareDescription();
    virtual void initHomepageRevisionText();
    virtual void initDevelopersList();
    virtual void initMainDevelopersList();
    virtual void initLicenseText();
protected:
    VoreenQtMainWindow* mainWindow_;        ///< pointer to the main window
    QString imagePath_;                     ///< current image path
    bool initialized_;                      ///< stores if the box has been initialized

    std::string softwareDescription_;       ///< stores the shown string
    std::string homepageRevisionText_;      ///< stores the shown string
    QStringList developers_;                ///< stores the shown string
    QStringList mainDevelopers_;            ///< stores the shown string
    std::string licenseString_;             ///< stores the shown string
    //---------------------
    //  help function
    //---------------------
protected:
    /**
     * Converts the StringList into a comma separated string.
     */
    std::string convertStringListToString(const QStringList& list);
};

} // namespace

#endif // VRN_ABOUTBOXBASE_H
