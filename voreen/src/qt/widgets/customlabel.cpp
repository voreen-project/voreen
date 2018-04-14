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

#include "voreen/qt/widgets/customlabel.h"
#include "voreen/qt/widgets/property/qpropertywidget.h"

//#include "voreen/core/properties/transfuncproperty.h"

#include <QContextMenuEvent>
#include <QFont>
#include <QLineEdit>
#include <QMenu>
#include <QWidget>

namespace voreen {

CustomLabel::CustomLabel(const char* text, QPropertyWidget* pw, QWidget* parent,
                         Qt::WindowFlags f, bool disallowWordwrap, bool editable)
    : QLabel(text, parent, f)
    , propertyWidget_(pw)
    , disallowWordwrap_(disallowWordwrap)
    , editable_(editable)
{
    init();
}

CustomLabel::CustomLabel(QPropertyWidget* pw, QWidget* parent, Qt::WindowFlags f, bool editable)
    : QLabel(parent, f)
    , propertyWidget_(pw)
    , disallowWordwrap_(false)
    , editable_(editable)
{
    init();
}

CustomLabel::~CustomLabel() {
    // disconnect from widget
    if (propertyWidget_)
        propertyWidget_->nameLabel_ = 0;
}

void CustomLabel::init() {
    initFont();
    processText();
    edit_ = new QLineEdit(this);
    edit_->hide();
    connect(edit_, SIGNAL(editingFinished()), this, SLOT(editingFinished()));
}

void CustomLabel::editingFinished() {
    setText(edit_->text());

    if(propertyWidget_ != 0) {
        propertyWidget_->setPropertyGuiName(edit_->text().toStdString());
    }
    edit_->hide();
}

void CustomLabel::contextMenuEvent(QContextMenuEvent* e) {
    if (editable_) {
        QMenu* men = new QMenu(this);
        men->addAction("Rename");
        men->addAction("Set Default");
        QAction* ac = men->exec(e->globalPos());
        if(ac != 0){
            if(ac->iconText().compare("Rename") == 0){
                edit_->setText(text());
                edit_->setFocus();
                edit_->setCursorPosition(edit_->text().length());
                edit_->resize(size());
                edit_->show();
            } else {
                propertyWidget_->getProperty()->reset();
                //if(dynamic_cast<TransFuncProperty*>(propertyWidget_->getProperty()))
                //    propertyWidget_->getProperty()->invalidate(); //y?
            }
        }
    }
}

void CustomLabel::initFont() {
    QFontInfo fontInfo(font());
    setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
    setToolTip(text());
}
void CustomLabel::processText() {
    //this should somehow reflect the size of a possible parent widget and therefore react on resize events
    if (text().length() > 20 && !disallowWordwrap_)
        setWordWrap(true);
}

} // namespace
