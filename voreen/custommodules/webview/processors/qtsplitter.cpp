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

#include "qtsplitter.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/processors/processorwidget.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

namespace voreen {

const std::string QtSplitter::loggerCat_("voreen.QtSplitter");

QtSplitter::QtSplitter()
    : Processor()
    , orientation_("orientation", "Orientation")
    , widgets_("widgets", "Widgets", false)
    , updateWidgets_("updateWidgets", "Update widgets")
{
    addProperty(orientation_);
    orientation_.addOption("horizontal", "Horizontal", HORIZONTAL);
    orientation_.addOption("vertical", "Vertical", VERTICAL);
    ON_CHANGE_LAMBDA(orientation_, [this] { getProcessorWidget()->updateFromProcessor(); });

    addProperty(widgets_);
    ON_CHANGE_LAMBDA(widgets_, [this] { getProcessorWidget()->updateFromProcessor(); });

    addProperty(updateWidgets_);
    ON_CHANGE(updateWidgets_, QtSplitter, networkChanged);
}

QtSplitter::~QtSplitter() {
}

Processor* QtSplitter::create() const {
    return new QtSplitter();
}

void QtSplitter::initialize() {
    Processor::initialize();

    if(auto* app = VoreenApplication::app()) {
        auto* network = app->getNetworkEvaluator()->getProcessorNetwork();
        network->addObserver(this);
    }

    networkChanged();
}

void QtSplitter::deinitialize() {

    if(auto* app = VoreenApplication::app()) {
        auto* network = app->getNetworkEvaluator()->getProcessorNetwork();
        network->removeObserver(this);
    }

    Processor::deinitialize();
}

bool QtSplitter::isReady() const {
    return true;
}

void QtSplitter::process() {
}

void QtSplitter::serialize(Serializer& s) const {
    Processor::serialize(s);
}
void QtSplitter::deserialize(Deserializer& s) {
    Processor::deserialize(s);
}

void QtSplitter::networkChanged() {

    if(!getProcessorWidget()) {
        return;
    }

    auto* network = VoreenApplication::app()->getNetworkEvaluator(this)->getProcessorNetwork();

    auto processors = std::vector<std::string>();
    for (auto* processor : network->getProcessors()) {
        if(processor != this && processor->getProcessorWidget() != nullptr) {
            processors.push_back(processor->getID());
        }
    }

    // Bookkeeping of already added instances (item and position).
    std::vector<std::pair<std::string, int>> keptWidgets;
    for(const std::string& processor : processors) {
        const auto& instances = widgets_.getInstances();
        for(int i=0, num=static_cast<int>(instances.size()); i<num; i++) {
            if(widgets_.getItems()[instances[i].getItemId()] == processor) {
                keptWidgets.emplace_back(std::pair<std::string, int>(processor, i));
            }
        }
    }

    // Setting items will reset instances.
    widgets_.blockCallbacks(true);
    widgets_.setItems(processors);

    // So we re-add them.
    for(const auto& item : keptWidgets) {
        widgets_.addInstance(item.first, item.second);
    }
    widgets_.blockCallbacks(false);

    getProcessorWidget()->updateFromProcessor();
}

QtSplitter::Orientation QtSplitter::getOrientation() const {
    return orientation_.getValue();
}

std::vector<std::string> QtSplitter::getInstances() const {
    std::vector<std::string> instances;
    for(const auto& instance : widgets_.getInstances()) {
        instances.push_back(instance.getName());
    }
    return instances;
}


} // namespace voreen
