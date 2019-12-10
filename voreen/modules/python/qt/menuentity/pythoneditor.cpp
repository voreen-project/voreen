/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

// include Python stuff first
#include "modules/python/pythonmodule.h"
#include "pythoneditor.h"
#include "../pythonplugin.h"

#include "voreen/qt/widgets/voreentoolwindow.h"

namespace voreen {

const std::string PythonEditor::loggerCat_ = "voreen.Python.PythonEditor";

PythonEditor::PythonEditor()
    : VoreenQtMenuEntity()
    , pythonWidget_(nullptr)
{}

PythonEditor::~PythonEditor() {
    delete pythonWidget_;
}

QWidget* PythonEditor::createWidget() const {
    pythonWidget_ = new PythonPlugin();
    return pythonWidget_;
}

void PythonEditor::initialize() {
    VoreenQtMenuEntity::initialize();

    if (PythonModule::getInstance()) {
        pythonWidget_->newScript();
    }
    else
        LWARNING("Python module class not instantiated");
}

void PythonEditor::deinitialize() {
    pythonWidget_->clearScript();

    VoreenQtMenuEntity::deinitialize();
}

} // namespace voreen
