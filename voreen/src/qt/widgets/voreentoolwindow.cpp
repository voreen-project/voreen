/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/voreentoolwindow.h"

#include <QAction>
#include <QApplication>
#include <QBoxLayout>
#include <QCloseEvent>
#include <QFrame>
#include <QPainter>
#include <iostream>

namespace voreen {

VoreenToolWindowTitle::VoreenToolWindowTitle(QDockWidget *parent, bool dockable)
    : QWidget(parent),
      dockable_(dockable)
{
    closeButton_ = QPixmap(":/qt/icons/closebutton.png");
    maximizeButton_ = QPixmap(":/qt/icons/maximizebutton.png");
    if (dockable_)
        undockButton_ = QPixmap(":/qt/icons/undockbutton.png");
}

QSize VoreenToolWindowTitle::sizeHint() const {
    return QSize(40,20);
}

QSize VoreenToolWindowTitle::minimumSizeHint() const {
    return sizeHint();
}

void VoreenToolWindowTitle::paintEvent(QPaintEvent* /*e*/) {
    QRect button_rect = rect();
    QPainter painter(this);

    QColor button_color = QColor(50, 50, 50);
    QColor m_shadow = QColor(Qt::black);
    int roundness = 0;

    // outline
    painter.setPen(QPen(QBrush(Qt::black), 2.0));
    QPainterPath outline;
    outline.addRoundRect(0, 0, button_rect.width(), button_rect.height(), roundness, roundness);
    painter.setOpacity(1.0);
    painter.drawPath(outline);

    // gradient
    QLinearGradient gradient(0, 0, 0, button_rect.height());
    gradient.setSpread(QGradient::ReflectSpread);
    gradient.setColorAt(0.0, button_color);
    gradient.setColorAt(0.4, m_shadow);
    gradient.setColorAt(0.6, m_shadow);
    gradient.setColorAt(1.0, button_color);

    QBrush brush(gradient);
    painter.setBrush(brush);
    painter.setPen(QPen(QBrush(button_color), 2.0));

    // main button
    QPainterPath painter_path;
    painter_path.addRoundRect(1, 1, button_rect.width() - 2, button_rect.height() - 2, roundness, roundness);
    painter.setClipPath(painter_path);

    painter.setOpacity(1.0);
    painter.drawRoundRect(1, 1, button_rect.width() - 2, button_rect.height() - 2, roundness, roundness);

    // glass highlight
    painter.setBrush(QBrush(Qt::white));
    painter.setPen(QPen(QBrush(Qt::white), 0.01));
    painter.setOpacity(0.30);
    painter.drawRect(1, 1, button_rect.width() - 2, (button_rect.height() / 2) - 1);

    // text
    QDockWidget* dockWidget = qobject_cast<QDockWidget*>(parentWidget());
    QString text = dockWidget->windowTitle();
    if (!text.isEmpty()) {
        painter.setFont(font());
        painter.setPen(Qt::white);
        painter.setOpacity(1.0);

        QRect text_rect = button_rect;
        text_rect.setLeft(8);
        text_rect.translate(0, 1);
        painter.drawText(text_rect, Qt::AlignLeft | Qt::AlignVCenter, text);
    }

    // draw logos
    QRect logoPos = button_rect;
    if (dockable_) {
        logoPos.setLeft(button_rect.right() - 32);
        logoPos.setTop(button_rect.top() + 4);
        painter.drawPixmap(logoPos.topLeft(), undockButton_);
        logoPos = button_rect;
    }

    logoPos.setLeft(button_rect.right() - 16);
    logoPos.setTop(button_rect.top() + 4);
    painter.drawPixmap(logoPos.topLeft(), closeButton_);

    QDockWidget *dw = qobject_cast<QDockWidget*>(parentWidget());
    if ( dw->isFloating() && (dw->maximumSize() != dw->minimumSize()) ) {
        logoPos.setLeft(button_rect.right() - 48);
        logoPos.setTop(button_rect.top() + 4);
        painter.drawPixmap(logoPos.topLeft(), maximizeButton_);
    }
}

void VoreenToolWindowTitle::mousePressEvent(QMouseEvent *event) {
    QPoint pos = event->pos();
    QRect rect = this->rect();

    QRect buttonRectClose = QRect(rect.right()-16, rect.top()+4, 14, 14);
    QRect buttonRectUndock = QRect(rect.right()-32, rect.top()+4, 14, 14);
    QRect buttonRectMaximize = QRect(rect.right()-48, rect.top()+4, 14, 14);
    if (buttonRectClose.contains(pos)) {
        event->accept();
        QDockWidget *dw = qobject_cast<QDockWidget*>(parentWidget());
        dw->close();
    } else if (buttonRectUndock.contains(pos) && dockable_) {
        event->accept();
        QDockWidget *dw = qobject_cast<QDockWidget*>(parentWidget());
        dw->setFloating(!dw->isFloating());
    } else if (buttonRectMaximize.contains(pos)) {
        QDockWidget *dw = qobject_cast<QDockWidget*>(parentWidget());
        if(!dw->isFloating() || (dw->maximumSize() == dw->minimumSize())) {
            event->ignore();
            return;
        }
        event->accept();
        if(dw->isMaximized())
            dw->showNormal();
        else
            dw->showMaximized();
    } else
        event->ignore();
}

void VoreenToolWindowTitle::mouseDoubleClickEvent(QMouseEvent *event) {
    // avoid docking via double click if desired
    if (!dockable_)
        event->accept();
    else
        event->ignore();
}

//----------------------------------------------------------------------------------------------------------------

VoreenToolWindow::VoreenToolWindow(QAction* action, QWidget* parent, QWidget* child, const QString& name, bool dockable)
    : QDockWidget(name, parent)
    , action_(action)
    , child_(child)
{
    if (action) {
        QString text = action->text();
        text.replace("&", ""); // remove menu accels
        setWindowTitle(text);
        setWindowIcon(action->icon());
    } else {
        setWindowTitle("");
    }
    setObjectName(name);

    QFrame* frame = new QFrame();
    frame->setFrameStyle(QFrame::StyledPanel);
    frame->setContentsMargins(1, 1, 1, 1);
    QVBoxLayout* vBox = new QVBoxLayout();
    vBox->addWidget(child);
    vBox->setContentsMargins(0, 0, 0, 0);
    frame->setLayout(vBox);
    setWidget(frame);

    // if we have a stylesheet we want the fancy window title, please
    if (!qApp->styleSheet().isEmpty())
        setTitleBarWidget(new VoreenToolWindowTitle(this, dockable));

    adjustSize();
    if (!dockable) {
        setAllowedAreas(Qt::NoDockWidgetArea);
        setFloating(true);
    }
    setVisible(false);

    // Connect action and widget visibility:
    // It is important to use triggered() instead of toggled() here, or else the widget will be
    // hidden when switching to a differen virtual desktop and back.
    if(action) {
        connect(action, SIGNAL(triggered(bool)), this, SLOT(setVisible(bool)));
        connect(toggleViewAction(), SIGNAL(toggled(bool)), action, SLOT(setChecked(bool)));
    }
}

void VoreenToolWindow::updateTitle(QString title, QIcon icon) {
    setWindowTitle(title.replace("&",""));
    setWindowIcon(icon);
}

QAction* VoreenToolWindow::action() const {
    return action_;
}

QWidget* VoreenToolWindow::child() const {
    return child_;
}

} // namespace
