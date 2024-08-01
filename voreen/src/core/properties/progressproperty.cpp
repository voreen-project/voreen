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

#include "voreen/core/properties/progressproperty.h"

namespace voreen {

ProgressProperty::ProgressProperty(const std::string& id, const std::string& guiText)
    : Property(id, guiText, Processor::VALID, Property::LOD_DEFAULT)
    , ProgressReporter()
    , message_("")
{}

ProgressProperty::ProgressProperty()
    : Property("", "", Processor::VALID, Property::LOD_DEFAULT)
    , ProgressReporter()
    , message_("")
{}

ProgressProperty::~ProgressProperty() {
}

Property* ProgressProperty::create() const {
    return new ProgressProperty();
}

void ProgressProperty::reset() {
    setProgress(0.f);
}

void ProgressProperty::setProgressMessage(const std::string& message) {
    message_ = message;
    updateWidgets();
}

std::string ProgressProperty::getProgressMessage() const {
    return message_;
}

void ProgressProperty::update() {
    updateWidgets();
}

}   // namespace
