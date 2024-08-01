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

#include "webviewwidget.h"
#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QMainWindow>
#include <QLabel>
#include <QCheckBox>
#include <QUrl>


namespace voreen {

const std::string WebViewWidget::loggerCat_("voreen.WebViewWidget");

WebViewWidget::WebViewWidget(QWidget* parent, WebViewProcessor* processor)
    : QProcessorWidget(processor, parent)
    , processor_(processor)
{
    tgtAssert(processor, "No WebView processor");

    setWindowTitle(QString::fromStdString(processor->getID()));
    resize(800, 480);
}

void WebViewWidget::initialize() {
    QProcessorWidget::initialize();

    webView_ = new QWebEngineView();

    QGridLayout* layout = new QGridLayout();
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(webView_);
    setLayout(layout);

    if(isVisible())
        show();

    initialized_ = true;
}

void WebViewWidget::updateFromProcessor() {
    webView_->load(QString::fromStdString(processor_->getURL()));
}

} //namespace voreen

