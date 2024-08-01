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

#ifndef VRN_PROGRESSBARGRAPHICSITEM_H
#define VRN_PROGRESSBARGRAPHICSITEM_H

#include "../nwebasegraphicsitem.h"
#include "voreen/core/io/progressbar.h"

#include "../../editor_settings.h"

#include <QTime>

namespace voreen {

class ProgressBarGraphicsItem : public NWEBaseGraphicsItem, public ProgressBar {
public:
    //constructor and destructor
    ProgressBarGraphicsItem(PortOwnerGraphicsItem* parent, const QPointF& center, qreal width, qreal height);
    virtual ~ProgressBarGraphicsItem();

     //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    int type() const {return UserTypesProgressBarGraphicsItem;}
    void updateNWELayerAndCursor() {};

    //style functions
    QRectF boundingRect() const;
    QPainterPath shape() const;
protected:
    void initializePaintSettings();
    void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);

public:
    void layoutChildItems(){};
protected:
    void createChildItems(){};
    void deleteChildItems(){};

    //---------------------------------------------------------------------------------
    //      progressbar functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    qreal getWidth() const;
    qreal getHeight() const;
    QColor getBackgroundColor1() const {return backgroundColor1_;}
    QColor getBackgroundColor2() const {return backgroundColor2_;}
    QColor getUpperForegroundColor1() const {return upperForegroundColor1_;}
    QColor getUpperForegroundColor2() const {return upperForegroundColor2_;}
    QColor getLowerForegroundColor() const {return lowerForegroundColor_;}


    void show();
    void hide();
    void forceUpdate();
    void update();

    void resize(const QPointF& center, qreal width, qreal height);

protected:
    qreal width_;
    qreal height_;
    QTime time_;

    QColor backgroundColor1_;
    QColor backgroundColor2_;
    QColor upperForegroundColor1_;
    QColor upperForegroundColor2_;
    QColor lowerForegroundColor_;
};

} // namespace

#endif // VRN_PROGRESSBARGRAPHICSITEM_H
