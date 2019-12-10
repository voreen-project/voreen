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

#ifndef VRN_VOLUMEURLISTPROPERTYWIDGET_H
#define VRN_VOLUMEURLISTPROPERTYWIDGET_H


#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/datastructures/volume/volumeram.h"

#include "voreen/qt/widgets/property/qpropertywidget.h"
#include "voreen/qt/widgets/volumeiohelper.h"

#include <QDialog>
#include <QCheckBox>
#include <QTreeWidget>

namespace voreen {


class VolumeURLListPropertyWidget : public QPropertyWidget {
    Q_OBJECT

public:
    VolumeURLListPropertyWidget(VolumeURLListProperty* volumeCollectionProp, QWidget* parent);

    /// Returns the null pointer, since this widget does not need a separate label.
    virtual CustomLabel* getOrCreateNameLabel() const;

protected slots:
    virtual void updateFromPropertySlot();
    virtual void changeFileWatchOption(bool checked);

private:
    void updateSelection();

    VolumeURLListProperty* urlListProperty_;
    VolumeIOHelper volumeIOHelper_;

    QTreeWidget* volumeTreeWidget_;
    QPushButton* loadButton_;
    QPushButton* clearButton_;
    QCheckBox* fileWatchCheckBox_;
    QCheckBox* selectAll_;

    static const std::string loggerCat_;

private slots:
    void volumeLoaded(const VolumeBase* handle);
    void clearVolumes();

    void itemSelected(QTreeWidgetItem*, int);
    void selectAll(bool);
};

} //namespace

#endif // VRN_VOLUMEURLISTPROPERTYWIDGET_H
