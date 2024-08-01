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

//header file
#include "cmlazylinker.h"
#include <limits>
#include <cstddef>

//we are in namespace voreen
namespace voreen {

CMLazyLinker::CMLazyLinker()
    : Processor()
    , intIn_("intIn", "Input: int", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , intOut_("intOut", "Output: int", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , floatIn_("floatIn", "Input: float", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , floatOut_("floatOut", "Output: float", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<float>::STATIC, Property::LOD_DEBUG)
    , cameraIn_("cameraIn", "Input: camera")
    , cameraOut_("cameraOut", "Output: camera")
    , intEvent_("intEvent", "Event trigger: int", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), Processor::VALID, NumericProperty<int>::STATIC, Property::LOD_DEBUG)

{
    addProperty(intIn_);
    addProperty(intOut_);
    addProperty(floatIn_);
    addProperty(floatOut_);
    addProperty(cameraIn_);
    addProperty(cameraOut_);
    addProperty(intEvent_);

    ON_CHANGE(intEvent_, CMLazyLinker, eventTriggered);
}

CMLazyLinker::~CMLazyLinker(){
}

void CMLazyLinker::initialize() {
    // call superclass function first
    Processor::initialize();
}

void CMLazyLinker::deinitialize() {
    Processor::deinitialize();
}
void CMLazyLinker::eventTriggered() {
    cameraOut_.set(cameraIn_.get());
    intOut_.set(intIn_.get());
    floatOut_.set(floatIn_.get());
}

void CMLazyLinker::process() {
}

} // namespace
