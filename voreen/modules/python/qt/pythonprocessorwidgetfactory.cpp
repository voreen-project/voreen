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

#include "pythonprocessorwidgetfactory.h"

#include "../processors/dynamicpythonprocessor.h"
#include "dynamicpythonwidget.h"

#include "voreen/qt/voreenapplicationqt.h"
#include "voreen/qt/mainwindow/voreenqtmainwindow.h"

namespace voreen {

ProcessorWidget* PythonProcessorWidgetFactory::createWidget(Processor* processor) const {

    if (!VoreenApplicationQt::qtApp()) {
        LERRORC("voreen.pythonProcessorWidgetFactory", "VoreenApplicationQt not instantiated");
        return 0;
    }
    QWidget* parent = VoreenApplicationQt::qtApp()->getMainWindow();

    if (dynamic_cast<DynamicPythonProcessor*>(processor))
        return new DynamicPythonWidget(parent, static_cast<DynamicPythonProcessor*>(processor));

    return 0;
}
} // namespace voreen
