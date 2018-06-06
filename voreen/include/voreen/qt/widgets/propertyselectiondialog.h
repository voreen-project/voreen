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

#ifndef VRN_CREATEPROPERTYSNAPSHOTDIALOG_H
#define VRN_CREATEPROPERTYSNAPSHOTDIALOG_H

#include <QDialog>

class QTreeWidget;

namespace voreen {

class ProcessorNetwork;
class Property;

class PropertySelectionDialog : public QDialog {
Q_OBJECT
public:
    PropertySelectionDialog(QList<Property*>& resultList, const ProcessorNetwork* network, QWidget* parent = 0);

public slots:
    void accept();

private:
    QList<Property*>& resultList_;
    const ProcessorNetwork* network_;
    QTreeWidget* treeWidget_;
};

} // namespace

#endif // VRN_CREATEPROPERTYSNAPSHOTDIALOG_H
