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

#include "voreen/qt/helpbrowser.h"

#include <QAction>
#include <QTextBrowser>
#include <QToolBar>

namespace voreen {

HelpBrowser::HelpBrowser(const QUrl& startPage, const QString& title, QWidget* parent)
    : QMainWindow(parent)
{
    browser_ = new QTextBrowser(this);
    browser_->setSource(startPage);
    browser_->createStandardContextMenu();
    setCentralWidget(browser_);

    setWindowIcon(QIcon(":/qt/icons/help.png"));

    setWindowTitle(title.isEmpty() ? tr("Voreen Help Browser") : title);

    createAndConnectToolbar();
}

void HelpBrowser::createAndConnectToolbar() {
    // create toolbar
    QToolBar* toolBar = new QToolBar(tr("Navigation"), this);

    // home
    QAction* homeAction = new QAction(QIcon(":/qt/icons/go-home.png"), tr("&Home"), this);
    homeAction->setStatusTip(tr("Show index"));
    homeAction->setToolTip(tr("Show first page"));
    homeAction->setShortcut(tr("Ctrl+H"));
    toolBar->addAction(homeAction);

    // backward
    QAction* backwardAction = new QAction(QIcon(":/qt/icons/go-previous.png"), tr("&Backward"), this);
    backwardAction->setStatusTip(tr("Show previous page"));
    backwardAction->setToolTip(tr("Show previous page"));
    backwardAction->setShortcut(tr("Ctrl+B"));
    toolBar->addAction(backwardAction);

    // forward
    QAction* forwardAction = new QAction(QIcon(":/qt/icons/go-next.png"), tr("&Forward"), this);
    forwardAction->setStatusTip(tr("Show next page"));
    forwardAction->setToolTip(tr("Show next page"));
    forwardAction->setShortcut(tr("Ctrl+F"));
    toolBar->addAction(forwardAction);

    // separator
    toolBar->addSeparator();

    // quit browser
    QAction* quitBrowserAction = new QAction(QIcon(":/qt/icons/exit.png"), tr("&Quit"), this);
    quitBrowserAction->setStatusTip(tr("Exit help"));
    quitBrowserAction->setToolTip(tr("Exit help"));
    quitBrowserAction->setShortcut(tr("Ctrl+W"));
    toolBar->addAction(quitBrowserAction);

    //setMenuBar(menuBar_);

    addToolBar(toolBar);

    // CONNECT ACTIONS
    // file menu
    connect(backwardAction, SIGNAL(triggered()), browser_, SLOT(backward()));
    connect(homeAction, SIGNAL(triggered()), browser_, SLOT(home()));
    connect(forwardAction, SIGNAL(triggered()), browser_, SLOT(forward()));
    connect(quitBrowserAction, SIGNAL(triggered()), this, SLOT(close()));
}

} // namespace voreen
