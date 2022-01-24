/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "webviewprocessor.h"

#include "voreen/core/processors/processorwidget.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

namespace voreen {

const std::string WebViewProcessor::loggerCat_("voreen.WebViewProcessor");

WebViewProcessor::WebViewProcessor()
    : Processor()
    , url_("url", "URL", "http://127.0.0.1:8050/")
{
    addProperty(url_);
    url_.setInstantUpdate(false);
    ON_CHANGE_LAMBDA(url_, [this] {
        if (getProcessorWidget()) {
            getProcessorWidget()->updateFromProcessor();
        }
    });
}

WebViewProcessor::~WebViewProcessor() {
}

Processor* WebViewProcessor::create() const {
    return new WebViewProcessor();
}

void WebViewProcessor::initialize() {
    Processor::initialize();
    getProcessorWidget()->updateFromProcessor();
}

void WebViewProcessor::deinitialize() {
    Processor::deinitialize();
}

bool WebViewProcessor::isReady() const {
    return true;
}

void WebViewProcessor::process() {
}

void WebViewProcessor::serialize(Serializer& s) const {
    Processor::serialize(s);
}
void WebViewProcessor::deserialize(Deserializer& s) {
    Processor::deserialize(s);
}


const std::string& WebViewProcessor::getURL() const {
    return url_.get();
}


} // namespace voreen
