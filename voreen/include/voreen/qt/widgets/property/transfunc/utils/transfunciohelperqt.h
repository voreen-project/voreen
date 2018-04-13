/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_TRANSFUNCIOHELPERQT_H
#define VRN_TRANSFUNCIOHELPERQT_H

#include "voreen/qt/voreenqtapi.h"

#include <QObject>
#include <QString>
#include <QStringList>

namespace voreen {

class TransFuncBase;

/**
 * Helper to encapsulate transfer function related IO functions.
 */
class VRN_QT_API TransFuncIOHelperQt : public QObject {
    Q_OBJECT
public:
    /**
     * Opens a filedialog and returns the choosen filename.
     *
     * @param filter filter with endings of supported file formats
     * @return the choosen filename
     */
    static const QString getLoadFileName(QStringList filters);

    /**
     * Opens a savefiledialog and returns the filename the user entered.
     *
     * @param filters filter with endings of supported file formats
     * @return filename the user entered.
     */
    static const QString getSaveFileName(QStringList filters);

    /**
    * Opens a Qt load dialog and loads the selected function into the param 'func'.
    * If an error occurs, the param will remain unchanged.
    * @param func function which will be replaced.
    * @return true, if the function has been changed.
    */
    static bool loadTransferFunction(TransFuncBase* func);

    /**
    * Opens a Qt load dialog and saves the the param 'func' to the selected file.
    * @param func function which will be saved.
    * @return true, if the function has been saved.
    */
    static bool saveTransferFunction(TransFuncBase* func);
};

} // namespace voreen

#endif // VRN_TRANSFUNCIOHELPERQT_H
