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

#ifndef VRN_EXPENDABLEHEADERBUTTON_H
#define VRN_EXPENDABLEHEADERBUTTON_H

#include <QFrame>

class QLabel;
class QToolButton;

namespace voreen {

class ExpandableHeaderButton : public QFrame {
Q_OBJECT
public:
    ExpandableHeaderButton(const QString& name, QWidget* parent = 0, bool expanded = false, bool userExpandanble = true);

    void setExpanded(bool expanded);
    void setUserExpandable(bool userExpandable);
    bool isExpanded() const;
    bool userExpandable() const;
    void updateNameLabel(const std::string& name);

    QFont fontStyle() const;
    void setFontStyle(QFont font);

signals:
    void toggled(bool expanded);

protected:
    void mousePressEvent(QMouseEvent* event);
    void updateState();

    QLabel* expandLabel_;
    QLabel* rendName_;
    bool expanded_;
    bool userExpandable_;
};

} // namespace

#endif // VRN_EXPENDABLEHEADERBUTTON_H
