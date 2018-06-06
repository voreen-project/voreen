/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_PROCESSORWIDGETFACTORY_H
#define VRN_PROCESSORWIDGETFACTORY_H

namespace voreen {

class Processor;
class ProcessorWidget;

/**
 * Interface for factories that create ProcessorWidgets.
 *
 * Usually, a module provides a ProcessorWidgetFactory for its processors.
 * However, since a processor does not necessarily need to be represented
 * by a specialized widget, providing a ProcessorWidgetFactory is optional.
 *
 * @see VoreenModule::registerProcessorWidgetFactory
 */
class ProcessorWidgetFactory {
public:
    virtual ~ProcessorWidgetFactory() {};
    virtual ProcessorWidget* createWidget(Processor* processor) const = 0;
};

} //namespace voreen

#endif
