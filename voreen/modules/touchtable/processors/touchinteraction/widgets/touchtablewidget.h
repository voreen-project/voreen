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

#ifndef VRN_TOUCHTABLEWIDGET_H
#define VRN_TOUCHTABLEWIDGET_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/genericcoprocessorport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "tgt/texturemanager.h"
#include  "voreen/core/voreenapplication.h"


namespace voreen {

    class TouchTableOverlay;

    /**
    * Abstract base class for widgets to be used with TouchTableOverlay (@see TouchTableOverlay).
    */
    class TouchTableWidget : public Processor {

    public:

        TouchTableWidget();

        virtual Processor* create() const = 0;

        virtual std::string getClassName() const;

        virtual std::string getCategory() const;

        virtual void initialize();

        virtual void deinitialize();

        virtual bool isReady() const;

        /**
        * Returns if the widget should be shown as active in the GUI.
        */
        virtual bool isActive() const;

        /**
        * Sets the active state of the widget.
        */
        virtual void setActive(bool active = true);

        /**
        * Called by TouchTableOverlay if a TouchPoint associated with this widget is released, ie. the widget has been pressed.
        */
        virtual void pressed() = 0;

        /**
        * Returns the current position of the center of the widget in canvas coordinates.
        */
        virtual tgt::vec2 getPosition() const;

        /**
        * Sets the widget center position in canvas coordinates.
        */
        virtual void setPosition(tgt::vec2 pos);

        /**
        * Returns the last valid position this widget was placed (in canvas coordinates).
        * Might be the same as the current position, if the current position is valid (ie. the widget is not currently moved or there is no collision).
        */
        virtual tgt::vec2 getOldPosition() const;

        /**
        * Sets the last valid position the widget was placed (in canvas coordinates).
        */
        virtual void setOldPosition(tgt::vec2 pos);

        /**
        * This returns the texture that should be rendered above the general widget layout by the TouchTableOverlay.
        */
        virtual tgt::Texture* getSymbolTexture() const;

        /**
        * True, if the widget is currently moved.
        */
        virtual bool isInMotion() const;

        /**
        * Set to true, if the widget is currently moved.
        */
        virtual void setInMotion(bool inMotion = true);

        /**
        * Returns true if the current position of the widget is invalid, ie. there is a collision.
        */
        virtual bool isPositionInvalid() const;

        /**
        * Set the current position of the widget invalid, e.g. if it is colliding with another widget while moving.
        */
        virtual void setPositionInvalid(bool invalid = true);

        /**
        * Reloads the image if file path has changed.
        */
        virtual void invalidate(int inv = INVALID_RESULT);

        /**
        * Set to true if the widget is currently placed on the canvas (ie. table), set to false if the widget is placed in the menu of the overlay.
        * @see TouchTableOverlay
        */
        virtual void setOnTable(bool t = true);

        /**
        * True, if the widget is currently placed on the table, false if it is placed in the menu.
        */
        virtual bool isOnTable() const;

        /**
         * Sets a pointer to the overlay this widget is connected to.
         */
        virtual void setOverlay(TouchTableOverlay* overlay);

    protected:

        /**
        * @see Processor
        */
        virtual void setDescriptions();

        /**
        * Check if coprocessor port has changed and sets the overlay to 0 if necessary.
        */
        virtual void process();

        GenericCoProcessorPort<TouchTableWidget> port_;

        BoolProperty isActive_;             ///< true if the widget should be rendered as activated by TouchTableOverlay
        BoolProperty inMotion_;             ///< true if the widget is currently moved
        BoolProperty invalidPosition_;      ///< true if the widget is currently in collision, ie. its position is invalid
        BoolProperty isOnTable_;            ///< true if the widget has been placed outside the menu on the table

        IntVec2Property position_;          ///< property for widget position in screen coordinates
        IntVec2Property oldPosition_;       ///< last valid position the widget has been

        FileDialogProperty symbolFile_;     ///< file dialog to select the image file used as a symbol texture

        tgt::Texture* symbolTex_;           ///< the symbol texture that is rendered above the general widget layout by the TouchTableOverlay

        TouchTableOverlay* overlay_;        ///< the overlay processor this widget os associated to (has to be set by overlay)

    };

} // namespace

#endif // VRN_TOUCHTABLEWIDGET_H
