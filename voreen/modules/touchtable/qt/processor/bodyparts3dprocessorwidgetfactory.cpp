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

#include "bodyparts3dprocessorwidgetfactory.h"

#include "bodyparts3dwidget.h"

#include "voreen/qt/voreenapplicationqt.h"
#include <QWidget>
#include <QMainWindow>

namespace voreen{

ProcessorWidget* BodyParts3DProcessorWidgetFactory::createWidget(Processor* processor) const{
    if(!VoreenApplicationQt::qtApp()){
        LERRORC("voreen.BodyParts3DProcessorWidgetFactory","VoreenApplicationQt not instantiated");
        return 0;
    }
    QWidget* parent = VoreenApplicationQt::qtApp()->getMainWindow();

    if(dynamic_cast<BodyParts3DSource*>(processor))
        return new BodyParts3DWidget(parent, static_cast<BodyParts3DSource*>(processor));

    return 0;
}

} // namespace voreen
