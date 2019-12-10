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

#ifndef VRN_PROCESSORWIDGET_H
#define VRN_PROCESSORWIDGET_H

#include "voreen/core/voreencoreapi.h"
#include <iostream>
#include "tgt/vector.h"

namespace voreen {

class Processor;

class VRN_CORE_API ProcessorWidget {
public:

    ProcessorWidget(Processor* processor);
    virtual ~ProcessorWidget();

    /**
     * Override this method for initializations that cannot or should not be done
     * in the constructor, e.g. OpenGL initializations. The function is called by
     * the owning processor's initialize() function.
     *
     * @note Call the superclass's function when overriding it,
     *       since it restores the window geometry.
     */
    virtual void initialize();

    /**
     * This function is called when the processor's state has changed.
     * Override it for updating the widget accordingly.
     */
    virtual void updateFromProcessor();

    /**
     * This function is called when the processor's name has changed.
     * Override it for updating the widget accordingly.
     */
    virtual void processorNameChanged();

    /**
     * Changes the visibility of the widget.
     * Has to be implemented by the derived class.
     */
    virtual void setVisible(bool) = 0;

    /**
     * Returns whether this widget is visible.
     * Has to be implemented by the derived class.
     */
    virtual bool isVisible() const = 0;

    /**
     * Resizes the widget.
     * Has to be implemented by the derived class.
     */
    virtual void setSize(int,int) = 0;

    /**
     * Returns the widget's size.
     * Has to be implemented by the derived class.
     */
    virtual tgt::ivec2 getSize() const = 0;

    /**
     * Moves the widget to the specified position.
     * Has to be implemented by the derived class.
     */
    virtual void setPosition(int,int) = 0;

    /**
     * Returns the widget's position.
     * Has to be implemented by the derived class.
     */
    virtual tgt::ivec2 getPosition() const = 0;

    /**
     * Saves the widget's geometry to processor's meta data for serialization.
     * Override for reacting to show events.
     */
    virtual void onShow();

    /**
     * Saves the widget's geometry to processor's meta data for serialization.
     * Override for reacting to hide events.
     */
    virtual void onHide();

    /**
     * Saves the widget's geometry to processor's meta data for serialization.
     * Override for reacting to move events.
     */
    virtual void onMove();

    /**
     * Saves the widget's geometry to processor's meta data for serialization.
     * Override for reacting to resize events.
     */
    virtual void onResize();

    /// Returns whether this widget has been initialized.
    bool isInitialized() const;

    /// Returns the porcessor
    const Processor* getProcessor() const;
protected:

    /// Saves the widget's geometry to the processor's meta data for serialization.
    virtual void saveGeometryToMeta();

    /// Restores the widget's geometry from the processor's meta data.
    virtual void restoreGeometryFromMeta();

    /// The processor the widget belongs to.
    Processor* processor_;

    /// Indicates whether this widget has been initialized.
    bool initialized_;

};

} //namespace voreen

#endif
