/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/expandableheaderbutton.h"

#include <QHBoxLayout>
#include <QLabel>
#include <QMouseEvent>
#include <QToolButton>

#include <tgt/assert.h>

namespace voreen {

ExpandableHeaderButton::ExpandableHeaderButton(const QString& name, QWidget* parent, bool expanded, bool expandable)
    : QFrame(parent)
    , expanded_(expanded)
    , userExpandable_(expandable)
{
    QHBoxLayout* frameLayout = new QHBoxLayout(this);
    frameLayout->setSpacing(0);
    frameLayout->setMargin(0);
    setFrameShape(QFrame::StyledPanel);

    expandLabel_ = new QLabel();
    expandLabel_->setMinimumWidth(20);
    expandLabel_->setMaximumWidth(20);
    frameLayout->addWidget(expandLabel_);

    // replace line breaks in proc name by spaces
    QString titleSingleLine = name;
    titleSingleLine.replace("\r\n", " ");
    titleSingleLine.replace("\r", " ");
    titleSingleLine.replace("\n", " ");

    rendName_ = new QLabel();
    rendName_->setAlignment(Qt::AlignCenter);
    rendName_->setText(titleSingleLine);
    frameLayout->addSpacing(2);
    frameLayout->addWidget(rendName_);

    updateState();
}

void ExpandableHeaderButton::mousePressEvent(QMouseEvent* event) {
    event->accept();
    if (userExpandable_ && (event->button() == Qt::LeftButton)) {
        setExpanded(!isExpanded());
        emit toggled(isExpanded());
    }
}

void ExpandableHeaderButton::setExpanded(bool expanded) {
    expanded_ = expanded;
    updateState();
}

bool ExpandableHeaderButton::isExpanded() const {
    return expanded_;
}

void ExpandableHeaderButton::updateState() {
    if (expanded_)
        expandLabel_->setPixmap(QPixmap(":/qt/icons/expand-minus.png"));
    else
        expandLabel_->setPixmap(QPixmap(":/qt/icons/expand-plus.png"));

    expandLabel_->setVisible(userExpandable_);
}

void ExpandableHeaderButton::setUserExpandable(bool userExpandable) {
    userExpandable_ = userExpandable;
    updateState();
}

bool ExpandableHeaderButton::userExpandable() const {
    return userExpandable_;
}

void ExpandableHeaderButton::updateNameLabel(const std::string& name) {
    QString titleSingleLine = QString::fromStdString(name);
    titleSingleLine.replace("\r\n", " ");
    titleSingleLine.replace("\r", " ");
    titleSingleLine.replace("\n", " ");
    rendName_->setText(titleSingleLine);
}

QFont ExpandableHeaderButton::fontStyle() const {
    tgtAssert(rendName_, "no name label");
    return rendName_->font();
}

void ExpandableHeaderButton::setFontStyle(QFont font) {
    tgtAssert(rendName_, "no name label");
    rendName_->setFont(font);
}


} // namespace
