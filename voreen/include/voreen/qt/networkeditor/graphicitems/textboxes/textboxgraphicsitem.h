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

#ifndef VRN_TEXTBOXGRAPHICSITEM_H
#define VRN_TEXTBOXGRAPHICSITEM_H

#include "textboxbasegraphicsitem.h"

namespace voreen {

    /**
     * This is the class used to represent documentation notes in the network editor.
     */
class TextBoxGraphicsItem : public TextBoxBaseGraphicsItem {
Q_OBJECT
public:
    //constructor + destructor
    TextBoxGraphicsItem(NetworkEditor* nwe);
    virtual ~TextBoxGraphicsItem();
    virtual int type() const {return UserTypesTextBoxGraphicsItem;}

    virtual void setContextMenuActions();
protected:
    virtual void initializePaintSettings();
    virtual void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
};

} // namespace

#endif // VRN_TEXTBOXGRAPHICSITEM_H




